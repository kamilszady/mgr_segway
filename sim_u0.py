"""
Symulacja UnTrans: model nieliniowy działa normalnie, ALE gdy |psi| przekroczy próg,
to BLOKUJEMY TYLKO:
  - psi (odchylenie)
  - psi_dot (pochodną odchylenia)
Pozostałe stany (theta, theta_dot, phi, phi_dot) ewoluują dalej wg modelu.

Wykorzystuje istniejącą klasę: UnTransNonlinearModel.derivatives()

Wymagania:
  pip install numpy scipy matplotlib
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt

# Import Twoich już zdefiniowanych klas (zmień nazwę modułu na swoją)
from untrans_model import UnTransNonlinearModel

psi_limit_deg: float = 102.0          # limit |psi| (próg "upadku"/ograniczenia)
post_limit_extra_sec: float = 2.0    # ile czasu po przekroczeniu limitu symulować dalej

@dataclass
class PsiClampOnly:
    psi_limit_deg: float = 102.0          # limit |psi| (próg "upadku"/ograniczenia)
    post_limit_extra_sec: float = 2.0    # ile czasu po przekroczeniu limitu symulować dalej


def simulate_clamp_psi_only(
    model: UnTransNonlinearModel,
    x0: np.ndarray,
    t_end: float,
    *,
    dt: float = 0.001,
    clamp: PsiClampOnly = PsiClampOnly(),
    u_override: Optional[Tuple[float, float]] = None,  # domyślnie u=0
):
    """
    RK4 + blokada tylko psi i psi_dot po przekroczeniu limitu.
    """
    x = np.asarray(x0, dtype=float).reshape(-1)
    if x.size != 6:
        raise ValueError("x0 must have length 6")

    psi_lim = float(clamp.psi_limit_deg) * np.pi / 180.0

    t = 0.0
    limited = False
    t_limit = None

    T = [t]
    Y = [x.copy()]

    def u_of(_t: float, _x: np.ndarray) -> Tuple[float, float]:
        return u_override if u_override is not None else (0.0, 0.0)

    def f(t_local: float, x_local: np.ndarray) -> np.ndarray:
        u_l, u_r = u_of(t_local, x_local)
        dx = np.asarray(model.derivatives(t_local, x_local, u_l=u_l, u_r=u_r), dtype=float)

        # Jeśli limit aktywny -> blokuj tylko psi i psi_dot
        if limited:
            dx[2] = 0.0  # psi_dot = 0 w sensie pochodnej psi
            dx[3] = 0.0  # psi_ddot = 0
        return dx

    def rk4_step(t_local: float, x_local: np.ndarray) -> np.ndarray:
        k1 = f(t_local, x_local)
        k2 = f(t_local + dt/2, x_local + dt*k1/2)
        k3 = f(t_local + dt/2, x_local + dt*k2/2)
        k4 = f(t_local + dt,   x_local + dt*k3)
        return x_local + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)

    # docelowy czas symulacji
    t_target = float(t_end)

    while t < t_target - 1e-12:
        x_next = rk4_step(t, x)

        # aktywacja limitu (raz)
        if (not limited) and abs(x_next[2]) >= psi_lim:
            limited = True
            t_limit = t + dt
            # "projekcja" tylko psi i psi_dot
            x_next[2] = np.sign(x_next[2]) * psi_lim if x_next[2] != 0 else psi_lim
            x_next[3] = 0.0
            # wydłuż symulację o extra czas po ograniczeniu
            t_target = max(t_target, t_limit + clamp.post_limit_extra_sec)

        # jeśli limit już aktywny, pilnuj żeby psi i psi_dot pozostały zablokowane
        if limited:
            x_next[2] = np.sign(x_next[2]) * psi_lim if x_next[2] != 0 else psi_lim
            x_next[3] = 0.0

        t += dt
        x = x_next
        T.append(t)
        Y.append(x.copy())

    T = np.asarray(T, dtype=float)
    Y = np.stack(Y, axis=1)  # (6,N)
    return T, Y, limited, t_limit


def plot_all_states(
    t: np.ndarray,
    y: np.ndarray,
    limited: bool,
    t_limit: Optional[float],
    psi_lim_deg: float,
):
    """
    Rysuje wszystkie 6 zmiennych stanu na jednej figurze (subplot pod subplotem).
    """
    labels = [
        r"$\theta$ [rad]",
        r"$\dot{\theta}$ [rad/s]",
        r"$\psi$ [rad]",
        r"$\dot{\psi}$ [rad/s]",
        r"$\phi$ [rad]",
        r"$\dot{\phi}$ [rad/s]",
    ]

    fig, axes = plt.subplots(6, 1, sharex=True, figsize=(10, 12))
    fig.suptitle(
        f"UnTrans – wszystkie zmienne stanu\n"
        f"Blokada tylko ψ, ψ̇ przy |ψ| = {psi_lim_deg}°",
        fontsize=12
    )

    for i, ax in enumerate(axes):
        ax.plot(t, y[i], linewidth=1.5)
        ax.set_ylabel(labels[i])
        ax.grid(True)

        if limited and t_limit is not None:
            ax.axvline(
                t_limit,
                linestyle="--",
                linewidth=1.2,
                color="r",
                label="limit ψ" if i == 0 else None
            )

        if i == 2:  # wykres psi — pokaż granice
            psi_lim = psi_lim_deg * np.pi / 180.0
            ax.axhline(+psi_lim, linestyle=":", color="k")
            ax.axhline(-psi_lim, linestyle=":", color="k")

    axes[-1].set_xlabel("t [s]")

    if limited:
        axes[0].legend(loc="upper right")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

if __name__ == "__main__":
    model = UnTransNonlinearModel()

    deg = np.pi / 180.0

    # Warunki początkowe tak, żeby dobiło do limitu bez sterowania
    x0 = np.array([
        0.0,          # theta
        0.0,          # theta_dot
        10.0 * deg,   # psi
        30.0 * deg,   # psi_dot
        0.0,          # phi
        0.0           # phi_dot
    ])

    clamp = PsiClampOnly(psi_limit_deg=psi_limit_deg, post_limit_extra_sec=post_limit_extra_sec)

    t, y, limited, t_lim = simulate_clamp_psi_only(
      model,
      x0=x0,
      t_end=3.0,
      dt=0.001,
      clamp=clamp,
    )

    plot_all_states(
        t=t,
        y=y,
        limited=limited,
        t_limit=t_lim,
        psi_lim_deg=clamp.psi_limit_deg,
    )
