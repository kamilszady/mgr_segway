from __future__ import annotations

from dataclasses import dataclass
from math import sin, cos
from typing import Tuple, Iterable, Optional


@dataclass
class UnTransParams:
    """
    Parametry modelu nieliniowego UnTrans.

    Jednostki zgodne z manualem:
      m, M [kg], R,W,L [m], J* [kg*m^2], Rm [ohm], Kt/Ke [Nm/A], [Vs/rad], fm [N*m*s?] (w manualu jako współczynnik tarcia)
    """
    # Masa i geometria
    m: float = 0.32          # masa jednego koła [kg]
    M: float = 5.41          # masa "nadwozia" [kg]
    R: float = 0.15 / 2.0    # promień koła [m] (w manualu: 2R=0.15)
    W: float = 0.40          # rozstaw kół (szerokość) [m]
    L: float = 0.102         # wysokość środka masy [m]

    # Bezwładności
    Jw: float = 0.0013       # bezwładność koła [kg*m^2] (zidentyfikowana)
    Jpsi: float = 0.104      # bezwładność osi przechyłu [kg*m^2] (zidentyfikowana)
    Jphi: float = 0.0484     # bezwładność osi obrotu w planie [kg*m^2]
    Jm: float = 0.00119      # bezwładność silnik+przekładnia (z uwzgl. przełożenia) [kg*m^2]

    # Silnik i tarcie
    Rm: float = 1.0          # rezystancja uzwojenia [ohm] (w manualu: RDC)
    Kt: float = 0.025        # stała momentu [Nm/A] (w manualu: Kt)
    Ke: float = 0.025        # stała SEM [Vs/rad] (w manualu: Ke)
    fm: float = 0.00024      # współczynnik tarcia (manual)

    # Grawitacja
    g: float = 9.81          # [m/s^2]


class UnTransNonlinearModel:
    """
    Nieliniowy model UnTrans:
      x = [theta, theta_dot, psi, psi_dot, phi, phi_dot]
      u = [u_l, u_r] gdzie u_* ~ PWM w zakresie [-1,1] (1 -> 100% PWM ~ 12V) wg manuala.

    Implementacja bazuje na równaniach (4.33–4.35) oraz (4.36–4.39).
    """

    def __init__(self, params: Optional[UnTransParams] = None):
        self.p = params or UnTransParams()

    # --- pomocnicze współczynniki z manuala ---
    def _d1_d2(self) -> Tuple[float, float]:
        # manual: d1 = Kl/Rm, d2 = Kl*Kb/Rm + fm
        # przyjęcie mapowania: Kl~Kt, Kb~Ke, Rm~RDC
        d1 = self.p.Kt / self.p.Rm
        d2 = (self.p.Kt * self.p.Ke) / self.p.Rm + self.p.fm
        return d1, d2

    def derivatives(self, t: float, x: Iterable[float], u_l: float, u_r: float) -> Tuple[float, float, float, float, float, float]:
        """
        Zwraca x_dot dla stanu x i sterowań (u_l, u_r).
        """
        theta, theta_dot, psi, psi_dot, phi, phi_dot = x

        p = self.p
        d1, d2 = self._d1_d2()

        # --- współczynniki z (4.31–4.32) ---
        c1 = ((2.0 * p.m + p.M) * (p.R ** 2) + 2.0 * p.Jw + 2.0 * p.Jm)
        c2 = (p.M * (p.L ** 2) + p.Jpsi + 2.0 * p.Jm)
        f1 = (p.M * p.R * p.L * cos(psi) - 2.0 * p.Jm)
        f2 = (p.M * p.R * p.L * (psi_dot ** 2) * sin(psi))
        f3 = (p.M * (p.L ** 2) * (phi_dot ** 2) * sin(psi) * cos(psi) + p.M * p.g * p.L * sin(psi))

        denom = (c1 * c2 - f1 * f1)
        if abs(denom) < 1e-9:
            # zabezpieczenie numeryczne: bardzo rzadko, ale może się pojawić przy ekstremalnych stanach
            denom = 1e-9 if denom >= 0 else -1e-9

        # --- siły uogólnione (4.36–4.38) ---
        # F_theta = d1(u_r + u_l) + 2 d2 (psi_dot - theta_dot)
        F_theta = d1 * (u_r + u_l) + 2.0 * d2 * (psi_dot - theta_dot)

        # F_psi = - d1(u_r + u_l) + 2 d2 (theta_dot - psi_dot)
        F_psi = -d1 * (u_r + u_l) + 2.0 * d2 * (theta_dot - psi_dot)

        # F_phi = - d1 * (W/(2R)) (u_r - u_l) - (W^2/(2R^2)) d2 * phi_dot
        F_phi = -d1 * (p.W / (2.0 * p.R)) * (u_r - u_l) - ((p.W ** 2) / (2.0 * (p.R ** 2))) * d2 * phi_dot

        # --- przyspieszenia (4.33–4.35) ---
        theta_ddot = (c2 * (F_theta + f2) - f1 * (F_psi + f3)) / denom
        psi_ddot   = (c1 * (F_psi + f3) - f1 * (F_theta + f2)) / denom

        denom_phi = (0.5 * p.m * (p.W ** 2) +
                     p.M * (p.L ** 2) * (sin(psi) ** 2) +
                     (p.W ** 2) / (2.0 * (p.R ** 2)) * (p.Jw + p.Jm) +
                     p.Jphi)

        if abs(denom_phi) < 1e-12:
            denom_phi = 1e-12

        phi_ddot = (F_phi - 2.0 * p.M * (p.L ** 2) * phi_dot * psi_dot * sin(psi) * cos(psi)) / denom_phi

        # --- pochodne stanu ---
        return (
            theta_dot,
            theta_ddot,
            psi_dot,
            psi_ddot,
            phi_dot,
            phi_ddot
        )

    def step_euler(self, x: Tuple[float, float, float, float, float, float], u_l: float, u_r: float, dt: float, t: float = 0.0):
        """
        Prosty krok całkowania Eulera (do szybkich testów).
        Do symulacji lepiej użyć solve_ivp (na dole przykład).
        """
        dx = self.derivatives(t, x, u_l, u_r)
        return tuple(xi + dt * dxi for xi, dxi in zip(x, dx))


# --- PRZYKŁAD UŻYCIA (symulacja) ---
if __name__ == "__main__":
    # Stan początkowy: lekko przechylony (psi)
    x0 = (0.0, 0.0, 0.05, 0.0, 0.0, 0.0)  # [theta, theta_dot, psi, psi_dot, phi, phi_dot]
    model = UnTransNonlinearModel()

    # Sterowanie (np. oba koła tak samo -> jazda do przodu/tyłu)
    u_l, u_r = 0.0, 0.0

    # Szybka symulacja Eulerem
    dt = 0.001
    t = 0.0
    x = x0
    for _ in range(2000):  # 2 sekundy
        x = model.step_euler(x, u_l, u_r, dt, t=t)
        t += dt

    print("Stan po 2s (Euler):", x)

    # Jeśli masz SciPy, lepiej tak:
    try:
        from scipy.integrate import solve_ivp

        def f(t, x):
            # tu możesz np. podać regulator:
            # u_l, u_r = ...
            return model.derivatives(t, x, u_l=0.0, u_r=0.0)

        sol = solve_ivp(f, t_span=(0.0, 2.0), y0=list(x0), max_step=0.005, rtol=1e-7, atol=1e-9)
        print("Stan po 2s (solve_ivp):", sol.y[:, -1])
    except ImportError:
        pass
