"""
UnTrans: stabilizacja + pozycja (servo) przez "PID-like" MIMO (state feedback + integratory),
nastawy z lokacji biegunów na modelu zlinearyzowanym w (x_op, u_op).

- Linearyzacja: różnice skończone -> A(6x6), B(6x2)
- Integratory: y = [theta, psi, phi] (możesz zmienić)
- Pole placement na układzie augmentowanym: [x; z] ma wymiar 6+ny
- Symulacja na nieliniowym modelu (RK4, stały krok dt)
- Wszystkie 6 zmiennych stanu na jednej figurze (6 wykresów pod sobą)

Wymagania:
  pip install numpy scipy matplotlib
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, Optional

import numpy as np
from scipy.signal import place_poles
import matplotlib.pyplot as plt

from untrans_model import UnTransNonlinearModel


# -----------------------------
# 1) Linearyzacja w punkcie pracy
# -----------------------------
def linearize_fd(
    model: UnTransNonlinearModel,
    x_op: np.ndarray,
    u_op: np.ndarray,
    *,
    eps_x: float = 1e-6,
    eps_u: float = 1e-6,
) -> Tuple[np.ndarray, np.ndarray]:
    x_op = np.asarray(x_op, dtype=float).reshape(6)
    u_op = np.asarray(u_op, dtype=float).reshape(2)

    def f(x: np.ndarray, u: np.ndarray) -> np.ndarray:
        return np.asarray(
            model.derivatives(0.0, x, u_l=float(u[0]), u_r=float(u[1])),
            dtype=float
        ).reshape(6)

    A = np.zeros((6, 6), dtype=float)
    B = np.zeros((6, 2), dtype=float)

    # A = df/dx
    for i in range(6):
        dx = np.zeros(6)
        dx[i] = eps_x
        fp = f(x_op + dx, u_op)
        fm = f(x_op - dx, u_op)
        A[:, i] = (fp - fm) / (2 * eps_x)

    # B = df/du
    for j in range(2):
        du = np.zeros(2)
        du[j] = eps_u
        fp = f(x_op, u_op + du)
        fm = f(x_op, u_op - du)
        B[:, j] = (fp - fm) / (2 * eps_u)

    return A, B


# --------------------------------------------
# 2) "PID-like" MIMO przez lokację biegunów
# --------------------------------------------
@dataclass
class PIDLikeDesign:
    """
    Definicja wyjść, które całkujemy (I-term).
    Domyślnie: theta, psi, phi => indeksy [0, 2, 4].
    """
    y_indices: Tuple[int, ...] = (0, 2, 4)   # theta, psi, phi
    poles: Tuple[complex, ...] = (
        # 6 biegunów dla x (balans szybszy, yaw/pozycja wolniejsze) + 3 bieguny dla integratorów
        -8, -9, -10,      # balans / szybkie tryby
        -2, -2.5, -3,     # wolniejsze (pozycja/yaw)
        -0.8, -1.0, -1.2  # integratory (wyraźnie wolniejsze!)
    )


def design_pidlike_poleplacement(
    A: np.ndarray,
    B: np.ndarray,
    design: PIDLikeDesign,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Projektuje regulator:
      u = -Kx x - Ki z
      z_dot = -C x   (dla x_ref=0; ogólnie z_dot = C(x_ref - x))

    Zwraca:
      Kx: (2,6)
      Ki: (2,ny)
      C:  (ny,6)
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)

    y_idx = list(design.y_indices)
    ny = len(y_idx)

    # C wybiera wyjścia do całkowania
    C = np.zeros((ny, 6), dtype=float)
    for k, i in enumerate(y_idx):
        C[k, i] = 1.0

    # Układ augmentowany:
    # x_dot = A x + B u
    # z_dot = C (x_ref - x) = -C x  (dla x_ref=0 w projekcie)
    #
    # [x_dot] = [A  0][x] + [B]u
    # [z_dot]   [-C 0][z]   [0]
    Aaug = np.block([
        [A,               np.zeros((6, ny))],
        [-C,              np.zeros((ny, ny))],
    ])
    Baug = np.vstack([B, np.zeros((ny, 2))])

    poles = design.poles
    if len(poles) != 6 + ny:
        raise ValueError(f"Potrzebujesz dokładnie {6+ny} biegunów, masz {len(poles)}")

    placed = place_poles(Aaug, Baug, poles)
    Kaug = placed.gain_matrix  # (2, 6+ny)

    Kx = Kaug[:, :6]
    Ki = Kaug[:, 6:]
    return Kx, Ki, C


# --------------------------------------------
# 3) Kontroler PID-like z anti-windup (MIMO)
# --------------------------------------------
class PIDLikeControllerMIMO:
    """
    u = u_op - Kx(x-x_ref) - Ki z
    z_dot = C(x_ref - x)

    Anti-windup: conditional integration dla z.
    Idea:
      - liczymy u_unsat
      - saturujemy -> u_sat
      - jeżeli saturacja, to dla każdego integratora i:
           jeśli kolumna Ki[:,i] "pcha" u dalej w stronę saturacji (w sensie znaku), zamrażamy z_i
    """

    def __init__(
        self,
        Kx: np.ndarray,
        Ki: np.ndarray,
        C: np.ndarray,
        *,
        u_min: float = -1.0,
        u_max: float = 1.0,
        z_limit: float = 50.0,
    ):
        self.Kx = np.asarray(Kx, dtype=float)
        self.Ki = np.asarray(Ki, dtype=float)
        self.C = np.asarray(C, dtype=float)

        self.u_min = float(u_min)
        self.u_max = float(u_max)
        self.z_limit = float(z_limit)

        self.ny = self.C.shape[0]
        self.z = np.zeros(self.ny, dtype=float)

    def reset(self):
        self.z[:] = 0.0

    def compute(
        self,
        x: np.ndarray,
        x_ref: np.ndarray,
        u_op: np.ndarray,
        dt: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Zwraca (u_sat, u_unsat) i aktualizuje z (z anti-windup).
        """
        x = np.asarray(x, dtype=float).reshape(6)
        x_ref = np.asarray(x_ref, dtype=float).reshape(6)
        u_op = np.asarray(u_op, dtype=float).reshape(2)

        e = x - x_ref  # stanowy błąd (dla Kx)
        y_err = self.C @ (x_ref - x)  # (ny,)

        # u przed saturacją
        u_unsat = u_op - (self.Kx @ e) - (self.Ki @ self.z)
        u_sat = np.clip(u_unsat, self.u_min, self.u_max)

        saturated = np.any(u_sat != u_unsat)

        # anti-windup dla integratorów:
        if saturated:
            # dla każdego z_i sprawdzamy czy integracja zwiększałaby saturację
            # (heurystyka: patrzymy na wpływ z_i na u = -Ki z, czyli du = -Ki[:,i] dz_i)
            dz = y_err * dt
            freeze = np.zeros(self.ny, dtype=bool)

            for i in range(self.ny):
                # przewidywana zmiana u od integracji z_i
                du_i = -self.Ki[:, i] * dz[i]  # (2,)
                # jeśli u jest na górze i du_i>0 -> pogarsza, jeśli u na dole i du_i<0 -> pogarsza
                worse_upper = (u_sat >= self.u_max - 1e-12) & (du_i > 0)
                worse_lower = (u_sat <= self.u_min + 1e-12) & (du_i < 0)
                if np.any(worse_upper | worse_lower):
                    freeze[i] = True

            dz[freeze] = 0.0
            self.z += dz
        else:
            self.z += y_err * dt

        self.z = np.clip(self.z, -self.z_limit, self.z_limit)

        # po aktualizacji z można przeliczyć u (opcjonalnie – daje spójność)
        u_unsat = u_op - (self.Kx @ e) - (self.Ki @ self.z)
        u_sat = np.clip(u_unsat, self.u_min, self.u_max)

        return u_sat.astype(float), u_unsat.astype(float)


# --------------------------------------------
# 4) Symulacja nieliniowa (RK4, stały dt)
# --------------------------------------------
@dataclass
class SimConfig:
    t_end: float = 30.0
    dt: float = 0.001


def rk4_step(model: UnTransNonlinearModel, t: float, x: np.ndarray, u: np.ndarray, dt: float) -> np.ndarray:
    def f(tt, xx):
        return np.asarray(model.derivatives(tt, xx, u_l=float(u[0]), u_r=float(u[1])), dtype=float)

    k1 = f(t, x)
    k2 = f(t + dt/2, x + dt*k1/2)
    k3 = f(t + dt/2, x + dt*k2/2)
    k4 = f(t + dt,   x + dt*k3)
    return x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)


def simulate_closed_loop_pidlike(
    model: UnTransNonlinearModel,
    ctrl: PIDLikeControllerMIMO,
    x0: np.ndarray,
    *,
    x_ref: np.ndarray,
    u_op: np.ndarray,
    cfg: SimConfig,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Zwraca:
      t: (N,)
      y: (6,N)
      u_log: (2,N)
    """
    x = np.asarray(x0, dtype=float).reshape(6)
    x_ref = np.asarray(x_ref, dtype=float).reshape(6)

    u_op = np.asarray(u_op, dtype=float).reshape(2)

    N = int(cfg.t_end / cfg.dt) + 1
    t_arr = np.linspace(0.0, cfg.t_end, N)

    y = np.zeros((6, N), dtype=float)
    u_log = np.zeros((2, N), dtype=float)

    ctrl.reset()

    for k, t in enumerate(t_arr):
        y[:, k] = x

        u, _ = ctrl.compute(x=x, x_ref=x_ref, u_op=u_op, dt=cfg.dt)
        u_log[:, k] = u

        if k < N - 1:
            x = rk4_step(model, t, x, u, cfg.dt)

    return t_arr, y, u_log


# --------------------------------------------
# 5) Wykres: wszystkie stany na jednej figurze
# --------------------------------------------
def plot_all_states_multi(
    trajectories: List[Tuple[np.ndarray, np.ndarray, str]],
    *,
    title: str,
):
    labels = [
        r"$\theta$ [rad]",
        r"$\dot{\theta}$ [rad/s]",
        r"$\psi$ [rad]",
        r"$\dot{\psi}$ [rad/s]",
        r"$\phi$ [rad]",
        r"$\dot{\phi}$ [rad/s]",
    ]

    fig, axes = plt.subplots(6, 1, sharex=True, figsize=(11, 13))
    fig.suptitle(title, fontsize=12)

    for i, ax in enumerate(axes):
        for (t, y, lab) in trajectories:
            ax.plot(t, y[i], linewidth=1.2, label=lab if i == 0 else None)
        ax.set_ylabel(labels[i])
        ax.grid(True)

    axes[-1].set_xlabel("t [s]")
    axes[0].legend(loc="upper right", fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


# =========================
# MAIN
# =========================
if __name__ == "__main__":
    model = UnTransNonlinearModel()

    # Punkt pracy (podajesz gdzie linearyzujesz)
    x_op = np.zeros(6, dtype=float)
    u_op = np.zeros(2, dtype=float)

    # Referencja: np. przejedź s_ref metrów -> theta_ref = s_ref / R
    s_ref = 0.5
    theta_ref = s_ref / model.p.R
    theta_ref = 0.0
    x_ref = np.array([theta_ref, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=float)

    # Linearyzacja
    A, B = linearize_fd(model, x_op=x_op, u_op=u_op)

    # Projekt regulatora PID-like przez lokację biegunów na układzie augmentowanym
    # design = PIDLikeDesign(
    #     y_indices=(0, 2, 4),  # całkujemy theta, psi, phi (możesz zmienić na np. (0,2) jeśli yaw Ci niepotrzebny)
    #     poles=(
    #         -10, -11, -12,     # szybkie (balans)
    #         -3, -3.5, -4,      # wolniejsze (jazda/yaw)
    #         -1.0, -1.2, -1.5   # integratory (wolniejsze od szybkich trybów!)
    #     )
    # )

    design = PIDLikeDesign(
        y_indices=(0,),        # tylko theta
        poles=(-10,-11,-12,-3,-3.5,-4,  -0.5)  # 6+1 biegunów
    )


    # y_indices=(0,4)
    # poles=(-10,-11,-12,-3,-3.5,-4,  -0.5,-0.7)


    Kx, Ki, C = design_pidlike_poleplacement(A, B, design)
    print("Kx=\n", Kx)
    print("Ki=\n", Ki)
    print("C=\n", C)

    ctrl = PIDLikeControllerMIMO(Kx, Ki, C, u_min=-1.0, u_max=1.0, z_limit=50.0)
    ctrl = PIDLikeControllerMIMO(Kx, Ki, C, u_min=-12.0, u_max=12.0, z_limit=50.0)

    # Zestaw warunków początkowych
    deg = np.pi / 180.0
    initials = [
        (np.array([0, 0,  2*deg, 0, 0, 0], float), "IC: psi=+2°"),
        (np.array([0, 0, -2*deg, 0, 0, 0], float), "IC: psi=-2°"),
        (np.array([0, 0,  4*deg, 0, 0, 0], float), "IC: psi=+4°"),
        (np.array([0, 0, -4*deg, 0, 0, 0], float), "IC: psi=-4°"),
        (np.array([0, 0,  4*deg, 8*deg, 0, 0], float), "IC: psi=+4°, psi_dot=+8°/s"),
        (np.array([0, 0, -4*deg,-8*deg, 0, 0], float), "IC: psi=-4°, psi_dot=-8°/s"),
    ]

    initials = [(np.array([0.0, 0.0, 5.5*deg, 0.0, 0.0, 0.0], float), "pionowanie (0)")]
    cfg = SimConfig()

    trajectories = []
    for x0, lab in initials:
        t, y, u_log = simulate_closed_loop_pidlike(
            model, ctrl, x0,
            x_ref=x_ref,
            u_op=u_op,
            cfg=cfg
        )
        trajectories.append((t, y, lab))

    plot_all_states_multi(
        trajectories,
        title="UnTrans – PID-like MIMO (Kx+Ki) z lokacji biegunów; linearyzacja w x*=0,u*=0"
    )
