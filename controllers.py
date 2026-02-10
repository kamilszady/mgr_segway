from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from scipy.linalg import solve_continuous_are


# =========================
#  Abstract base controller
# =========================

class Controller(ABC):
    """
    Abstrakcyjna klasa bazowa dla wszystkich kontrolerów.

    Konwencja:
      x = [theta, theta_dot, psi, psi_dot, phi, phi_dot]
      u = [u_l, u_r]
    """

    def __init__(self, u_min: float = -1.0, u_max: float = 1.0):
        self.u_min = float(u_min)
        self.u_max = float(u_max)

    @abstractmethod
    def reset(self) -> None:
        pass

    @abstractmethod
    def compute(self, t: float, x: np.ndarray, r: Optional[np.ndarray] = None) -> np.ndarray:
        pass

    def __call__(self, t: float, x: np.ndarray, r: Optional[np.ndarray] = None) -> np.ndarray:
        return self.saturate(self.compute(t, x, r))

    def saturate(self, u: np.ndarray) -> np.ndarray:
        return np.clip(u, self.u_min, self.u_max)

    @staticmethod
    def _vec(x: Sequence[float], n: int, name: str) -> np.ndarray:
        x = np.asarray(x, dtype=float).reshape(-1)
        if x.size != n:
            raise ValueError(f"{name} must have length {n}")
        return x


# ==========
#  PID
# ==========

@dataclass
class PIDGains:
    kp: float
    ki: float
    kd: float


class PIDController(Controller):
    """
    PID stabilizujący kąt przechyłu psi (x[2]) z ANTI-WINDUP (conditional integration / clamping).

    - Wyjście: u_l = u_r (sterowanie wspólne).
    - D-term: używa psi_dot (x[3]) -> dla stałego setpointu: d/dt(e) = -psi_dot.
    - Anti-windup:
        jeśli sterowanie jest w saturacji i błąd chce "pchać dalej" w tę samą stronę,
        to NIE aktualizujemy całki.

    Uwaga:
      - Zakładamy, że saturacja aktuatora to u_min..u_max (np. [-1,1] dla PWM).
      - integrator_limit ogranicza zakres całki jako dodatkowe zabezpieczenie.
    """

    def __init__(
        self,
        gains: PIDGains,
        dt: float,
        *,
        psi_idx: int = 2,
        psi_dot_idx: int = 3,
        integrator_limit: float = 10.0,
        u_min: float = -12.0,
        u_max: float = 12.0,
    ):
        super().__init__(u_min, u_max)

        self.g = gains
        self.dt = float(dt)
        if self.dt <= 0:
            raise ValueError("dt must be > 0")

        self.psi_idx = int(psi_idx)
        self.psi_dot_idx = int(psi_dot_idx)
        self.integrator_limit = float(integrator_limit)

        self._I = 0.0

    def reset(self) -> None:
        self._I = 0.0

    def compute(self, t: float, x: np.ndarray, r: Optional[np.ndarray] = None) -> np.ndarray:
        x = np.asarray(x, dtype=float).reshape(-1)

        psi = float(x[self.psi_idx])
        psi_dot = float(x[self.psi_dot_idx])

        psi_ref = 0.0 if r is None else float(np.asarray(r, dtype=float).reshape(-1)[self.psi_idx])
        e = psi_ref - psi

        # D: dla stałego ref: d/dt(e) = -psi_dot
        d = -psi_dot

        # 1) policz sterowanie "nienasycone" dla aktualnej całki
        u_unsat = (self.g.kp * e) + (self.g.ki * self._I) + (self.g.kd * d)

        # 2) saturacja aktuatora
        u_sat = float(np.clip(u_unsat, self.u_min, self.u_max))

        # 3) ANTI-WINDUP (conditional integration):
        # - jeśli jesteśmy w saturacji i błąd pcha dalej w stronę saturacji, to nie całkuj
        saturated = (u_sat != u_unsat)
        pushing_further = ((u_sat >= self.u_max and e > 0.0) or (u_sat <= self.u_min and e < 0.0))

        if not (saturated and pushing_further):
            self._I += e * self.dt
            self._I = float(np.clip(self._I, -self.integrator_limit, self.integrator_limit))

            # (opcjonalnie) przelicz u po aktualizacji całki, żeby mieć spójne wyjście
            u_unsat = (self.g.kp * e) + (self.g.ki * self._I) + (self.g.kd * d)
            u_sat = float(np.clip(u_unsat, self.u_min, self.u_max))

        return np.array([u_sat, u_sat], dtype=float)


# ==========
#  LQR (Q, R)
# ==========

class LQRController(Controller):
    """
    LQR: u = -K (x - r)

    K liczone automatycznie z macierzy:
      x_dot = A x + B u
      koszt: ∫ (xᵀ Q x + uᵀ R u) dt
    """

    def __init__(
        self,
        A: np.ndarray,
        B: np.ndarray,
        Q: np.ndarray,
        R: np.ndarray,
        *,
        x_ref: Optional[np.ndarray] = None,
        u_min: float = -1.0,
        u_max: float = 1.0,
    ):
        super().__init__(u_min, u_max)

        self.A = np.asarray(A, dtype=float)
        self.B = np.asarray(B, dtype=float)
        self.Q = np.asarray(Q, dtype=float)
        self.R = np.asarray(R, dtype=float)

        if self.A.shape != (6, 6):
            raise ValueError("A must be 6x6")
        if self.B.shape != (6, 2):
            raise ValueError("B must be 6x2")
        if self.Q.shape != (6, 6):
            raise ValueError("Q must be 6x6")
        if self.R.shape != (2, 2):
            raise ValueError("R must be 2x2")

        self.x_ref = np.zeros(6) if x_ref is None else self._vec(x_ref, 6, "x_ref")

        self.K = self._compute_lqr_gain()

    def _compute_lqr_gain(self) -> np.ndarray:
        """
        Rozwiązuje ciągłe równanie Riccatiego i liczy K.
        """
        P = solve_continuous_are(self.A, self.B, self.Q, self.R)
        K = np.linalg.inv(self.R) @ self.B.T @ P
        return K

    def reset(self) -> None:
        pass

    def compute(self, t: float, x: np.ndarray, r: Optional[np.ndarray] = None) -> np.ndarray:
        x = self._vec(x, 6, "x")
        r_vec = self.x_ref if r is None else self._vec(r, 6, "r")

        e = x - r_vec
        u = -self.K @ e
        return u
