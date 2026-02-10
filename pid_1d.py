"""
UnTrans: PID do stabilizacji (pionowanie) z nastawami wyznaczonymi przez lokację biegunów
na zlinearyzowanym modelu w zadanym punkcie pracy.

Założenia (zgodne z tym co już mamy):
- Masz klasę UnTransNonlinearModel z metodą:
    derivatives(t, x, u_l, u_r) -> (6,) pochodnych stanu
- Masz klasy kontrolerów, w tym PIDController (dziedziczy po Controller), który:
    * reguluje psi = x[2]
    * używa psi_dot = x[3] jako członu D (czyli klasyczny PID na ψ)
    * daje sterowanie wspólne: u_l = u_r = u
- Sterowanie u jest znormalizowanym PWM w [-1, 1]. :contentReference[oaicite:0]{index=0}

Co robi skrypt:
1) Linearyzuje nieliniowy model w punkcie (x_op, u_op) -> A,B (różnice skończone)
2) Sprowadza do SISO dla pionowania: wejście wspólne u = (u_l+u_r)/2, wyjście ψ
3) Buduje model 2-rzędu [ψ, ψ_dot] i dodaje całkę błędu (PID = P + D + I z całki)
4) Dobiera [k_p, k_d, k_i] przez place_poles na układzie augmentowanym
5) Symuluje nieliniowy model w pętli zamkniętej PID dla wielu warunków początkowych
6) Rysuje wszystkie 6 stanów na jednej figurze (6 wykresów pod sobą)

Wymagania:
  pip install numpy scipy matplotlib
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import place_poles
import matplotlib.pyplot as plt

# ===== UŻYWAMY TWOICH ISTNIEJĄCYCH KLAS =====
# Zmień nazwy modułów na swoje
from untrans_model import UnTransNonlinearModel
from controllers import PIDController, PIDGains


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


# ---------------------------------------------------------
# 2) Dobór PID przez lokację biegunów (na zlinearyzowanym)
# ---------------------------------------------------------
def pid_from_pole_placement_for_psi(
    A: np.ndarray,
    B: np.ndarray,
    *,
    desired_poles: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    """
    Projekt PID na ψ przy założeniu wejścia wspólnego:
        u = (u_l + u_r)/2
    i modelu 2-rzędu w stanie [ψ, ψ_dot].

    PID realizujemy jako:
        u = - (k_p * ψ + k_d * ψ_dot + k_i * z)
        z_dot = -ψ  (całka błędu, dla ψ_ref=0)

    Zatem dobieramy K = [k_p, k_d, k_i] przez place_poles na układzie augmentowanym (3x3).
    """

    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)

    # wejście wspólne (SISO): u = (u_l+u_r)/2
    Bc = 0.5 * (B[:, 0] + B[:, 1])  # (6,)

    # bierzemy tylko część dla ψ:
    # x_psi = [ψ, ψ_dot] = [x[2], x[3]]
    Apsi = A[np.ix_([2, 3], [2, 3])]  # 2x2
    Bpsi = Bc[[2, 3]].reshape(2, 1)   # 2x1

    # augmentacja o całkę:
    # x_aug = [ψ, ψ_dot, z]
    # z_dot = -ψ
    Aaug = np.array([
        [Apsi[0, 0], Apsi[0, 1], 0.0],
        [Apsi[1, 0], Apsi[1, 1], 0.0],
        [-1.0,       0.0,       0.0],
    ], dtype=float)


    # Wewnątrz pid_from_pole_placement_for_psi:
    Aaug = np.array([
        [Apsi[0, 0], Apsi[0, 1], 0.0],
        [Apsi[1, 0], Apsi[1, 1], 0.0],
        [1.0,        0.0,        0.0],  # Zmień z -1.0 na 1.0
    ], dtype=float)

    Baug = np.array([
        [Bpsi[0, 0]],
        [Bpsi[1, 0]],
        [0.0],
    ], dtype=float)

    # place poles: A_cl = Aaug - Baug*K
    placed = place_poles(Aaug, Baug, desired_poles)
    K = placed.gain_matrix.reshape(-1)  # [k_p, k_d, k_i]

    k_p, k_d, k_i = float(K[0]), float(K[1]), float(K[2])
    return k_p, k_d, k_i


# -----------------------------
# 3) Symulacja nieliniowa w pętli PID
# -----------------------------
@dataclass
class SimConfig:
    t_end: float = 100.0
    max_step: float = 0.002
    rtol: float = 1e-7
    atol: float = 1e-9


def simulate_pid_closed_loop(
    model: UnTransNonlinearModel,
    pid: PIDController,
    x0: np.ndarray,
    cfg: SimConfig,
) -> Tuple[np.ndarray, np.ndarray]:
    x0 = np.asarray(x0, dtype=float).reshape(6)

    # referencja pionu: psi_ref = 0
    x_ref = np.zeros(6, dtype=float)
    s_ref = 0.5  # odległość do przejechania (m)
    theta_ref = s_ref / model.p.R  # np. s_ref = 0.0
    theta_ref = 0.0
    x_ref = np.array([theta_ref, 0.0, 0.0, 0.0, 0.0, 0.0], float)

    def dyn(t: float, x: np.ndarray) -> np.ndarray:
        # PIDController -> u = [u_l, u_r], już z saturacją [-1,1]
        u = pid(t, x, x_ref)
        return np.asarray(model.derivatives(t, x, u_l=float(u[0]), u_r=float(u[1])), dtype=float)

    pid.reset()
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

    return sol.t, sol.y  # y: (6, N)


# -----------------------------
# 4) Wykres: wszystkie stany na jednej figurze
# -----------------------------
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

    # Punkt pracy do linearyzacji (pion, brak sterowania)
    x_op = np.zeros(6, dtype=float)
    u_op = np.zeros(2, dtype=float)

    # 1) linearyzacja
    A, B = linearize_fd(model, x_op=x_op, u_op=u_op, eps_x=1e-6, eps_u=1e-6)

    # 2) pole placement -> PID gains
    # Dobór biegunów: trzy stabilne, ujemne (im bardziej ujemne, tym szybciej, ale większa saturacja)
    desired_poles = (-2.0, -1.0, -1.5)
    desired_poles = (-1.0, -0.01, -0.001) # bardzo dobre, ale rózne znaki wzmocnień
    desired_poles = (-1.0, -0.1, -0.01)
    desired_poles = (-2.0, -2.1, -2.2) # bardzo dobre, wszystkie ujemne
    desired_poles = (-1.5, -1.6, -1.7) # spoko, wszystkie ujemne
    # desired_poles = (-0.5, -0.6, -0.7) # za wolne
    kp, kd, ki = pid_from_pole_placement_for_psi(A, B, desired_poles=desired_poles)

    print("PID from pole placement:")
    print(f"  desired poles: {desired_poles}")
    print(f"  kp = {kp:.6f}, ki = {ki:.6f}, kd = {kd:.6f}")

    # 3) Kontroler PID (Twoja klasa)
    # dt: używane przez integrator w PID; ustaw zgodnie z tym, jak chcesz "próbkować" integral (np. 0.001..0.01)
    dt_pid = 0.001
    pid = PIDController(
        gains=PIDGains(kp=kp, ki=ki, kd=kd),
        dt=dt_pid,
        psi_idx=2,
        psi_dot_idx=3,
        integrator_limit=50.0,
        u_min=-12.0,
        u_max=12.0
    )

    # 4) Zestaw warunków początkowych
    deg = np.pi / 180.0
    initials = [
        (np.array([0, 0,  2*deg, 0, 0, 0], float), "psi=+2°"),
        (np.array([0, 0, -2*deg, 0, 0, 0], float), "psi=-2°"),
        (np.array([0, 0,  4*deg, 0, 0, 0], float), "psi=+4°"),
        (np.array([0, 0, -4*deg, 0, 0, 0], float), "psi=-4°"),
        (np.array([0, 0,  4*deg, 8*deg, 0, 0], float), "psi=+4°, psi_dot=+8°/s"),
        (np.array([0, 0, -4*deg,-8*deg, 0, 0], float), "psi=-4°, psi_dot=-8°/s"),
    ]

    initials = [(np.array([0.0, 0.0, 0.5*deg, 0.0, 0.0, 0.0], float), "pionowanie (0)")]
    cfg = SimConfig()

    trajectories = []
    for x0, lab in initials:
        t, y = simulate_pid_closed_loop(model, pid, x0, cfg)
        trajectories.append((t, y, lab))

    plot_all_states_multi(
        trajectories,
        title="UnTrans – pionowanie PID; nastawy z lokacji biegunów (model zlinearyzowany w x*=0,u*=0)"
    )
