from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, Optional

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ====== UŻYWAMY TWOICH ISTNIEJĄCYCH KLAS ======
# Zmień nazwy modułów na swoje
from untrans_model import UnTransNonlinearModel
from controllers import LQRController


# -----------------------------
#  1) Linearyzacja w punkcie
# -----------------------------

def linearize_fd(
    model: UnTransNonlinearModel,
    x_op: np.ndarray,
    u_op: np.ndarray,
    *,
    eps_x: float = 1e-6,
    eps_u: float = 1e-6,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Numeryczna linearyzacja: A = df/dx, B = df/du w punkcie (x_op, u_op)
    dla modelu: x_dot = f(t, x, u_l, u_r)

    Zwraca:
      A: (6,6)
      B: (6,2)
    """
    x_op = np.asarray(x_op, dtype=float).reshape(6)
    u_op = np.asarray(u_op, dtype=float).reshape(2)

    def f(x: np.ndarray, u: np.ndarray) -> np.ndarray:
        return np.asarray(model.derivatives(0.0, x, u_l=float(u[0]), u_r=float(u[1])), dtype=float).reshape(6)

    f0 = f(x_op, u_op)

    # A = df/dx
    A = np.zeros((6, 6), dtype=float)
    for i in range(6):
        dx = np.zeros(6, dtype=float)
        dx[i] = eps_x
        fp = f(x_op + dx, u_op)
        fm = f(x_op - dx, u_op)
        A[:, i] = (fp - fm) / (2.0 * eps_x)

    # B = df/du
    B = np.zeros((6, 2), dtype=float)
    for j in range(2):
        du = np.zeros(2, dtype=float)
        du[j] = eps_u
        fp = f(x_op, u_op + du)
        fm = f(x_op, u_op - du)
        B[:, j] = (fp - fm) / (2.0 * eps_u)

    return A, B


# -------------------------------------------
#  2) Symulacja nieliniowa z kontrolerem LQR
# -------------------------------------------

@dataclass
class SimConfig:
    t_end: float = 20.0
    max_step: float = 0.002
    rtol: float = 1e-7
    atol: float = 1e-9


def simulate_closed_loop(
    model: UnTransNonlinearModel,
    controller: LQRController,
    x0: np.ndarray,
    *,
    x_ref: np.ndarray,
    u_op: np.ndarray,
    cfg: SimConfig,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Symuluje nieliniowy model w pętli zamkniętej:
      u = u_op + u_lqr, gdzie u_lqr = -K (x - x_ref)
    """
    x0 = np.asarray(x0, dtype=float).reshape(6)
    x_ref = np.asarray(x_ref, dtype=float).reshape(6)
    u_op = np.asarray(u_op, dtype=float).reshape(2)

    def dyn(t: float, x: np.ndarray) -> np.ndarray:
        # LQRController zwraca u = -K(x-r), z saturacją jeśli masz ustawione u_min/u_max
        u_lqr = controller(t, x, x_ref)  # shape (2,)
        u = u_op + u_lqr
        return np.asarray(model.derivatives(t, x, u_l=float(u[0]), u_r=float(u[1])), dtype=float)

    sol = solve_ivp(
        dyn,
        t_span=(0.0, cfg.t_end),
        y0=x0,
        method="RK45",
        max_step=cfg.max_step,
        rtol=cfg.rtol,
        atol=cfg.atol,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    return sol.t, sol.y  # y: (6,N)


# -------------------------------------------
#  3) Wykres: wszystkie stany na jednej figurze
# -------------------------------------------

def plot_all_states_multi(
    trajectories: List[Tuple[np.ndarray, np.ndarray, str]],
    *,
    title: str = "UnTrans – pionowanie LQR (nieliniowa symulacja, różne IC)",
):
    """
    trajectories: lista (t, y, label)
      t: (N,)
      y: (6,N)
    """
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

    # legenda tylko na pierwszym wykresie, żeby nie zalać figury
    axes[0].legend(loc="upper right", fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


# =========================
#  MAIN: konfiguracja eksperymentu
# =========================

if __name__ == "__main__":
    model = UnTransNonlinearModel()
    deg = np.pi / 180.0
    # --- Punkt pracy (x*, u*) do linearyzacji ---
    # Pionowanie: zwykle x* = 0, u* = 0
    x_op = np.zeros(6, dtype=float)
    u_op = np.zeros(2, dtype=float)

    # Referencja stabilizacji (pion)
    x_ref = np.zeros(6, dtype=float)
    s_ref = 0.5  # odległość do przejechania (m)
    theta_ref = s_ref / model.p.R  # np. s_ref = 0.0
    theta_ref = 0.0
    x_ref = np.array([theta_ref, 0.0, 0.0, 0.0, 0.0, 0.0], float)

    # --- Linearyzacja w podanym punkcie ---
    A, B = linearize_fd(model, x_op=x_op, u_op=u_op, eps_x=1e-6, eps_u=1e-6)

    # --- Dobór Q, R (przykład – dostosuj) ---
    # Mocno karz psi i psi_dot, umiarkowanie resztę
    Q = np.diag([0.5, 1.5, 20000.0, 0.1, 20.0, 750.0])
    R = np.diag([1000.0, 1000.0])

    # --- Kontroler LQR (Twoja klasa liczy K z A,B,Q,R) ---
    lqr = LQRController(A=A, B=B, Q=Q, R=R, x_ref=x_ref, u_min=-12.0, u_max=12.0) # dla +-12 lwpiej działa, logiczniej

    print("K (LQR gain) =\n", lqr.K)

    # --- Zestaw warunków początkowych (różne odchylenia) ---
    # initials = [
    #     (np.array([0, 0,  3*deg, 0, 0, 0], float), "psi=+3°"),
    #     (np.array([0, 0, -3*deg, 0, 0, 0], float), "psi=-3°"),
    #     (np.array([0, 0,  6*deg, 0, 0, 0], float), "psi=+6°"),
    #     (np.array([0, 0, -6*deg, 0, 0, 0], float), "psi=-6°"),
    #     (np.array([0, 0,  6*deg, 10*deg, 0, 0], float), "psi=+6°, psi_dot=+10°/s"),
    #     (np.array([0, 0, -6*deg, -10*deg, 0, 0], float), "psi=-6°, psi_dot=-10°/s"),
    # ]

    initials = [(np.array([0.0, 0.0, 5.0*deg, 0.0, 0.0, 0.0], float), "pionowanie (0)")]

    cfg = SimConfig()

    trajectories = []
    for x0, lab in initials:
        t, y = simulate_closed_loop(
            model=model,
            controller=lqr,
            x0=x0,
            x_ref=x_ref,
            u_op=u_op,
            cfg=cfg,
        )
        trajectories.append((t, y, lab))

    plot_all_states_multi(trajectories)
