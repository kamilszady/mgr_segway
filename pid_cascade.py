import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import place_poles
import matplotlib.pyplot as plt
from dataclasses import dataclass

# Zakładamy, że klasy są w tych samych plikach/folderach
from untrans_model import UnTransNonlinearModel
from controllers import PIDController, PIDGains

# -----------------------------
# 1) Linearyzacja (Twoja funkcja)
# -----------------------------
def linearize_fd(model, x_op, u_op, eps_x=1e-6, eps_u=1e-6):
    def f(x, u):
        return np.asarray(model.derivatives(0.0, x, u_l=float(u[0]), u_r=float(u[1]))).reshape(6)
    
    A, B = np.zeros((6, 6)), np.zeros((6, 2))
    for i in range(6):
        dx = np.zeros(6); dx[i] = eps_x
        A[:, i] = (f(x_op + dx, u_op) - f(x_op - dx, u_op)) / (2 * eps_x)
    for j in range(2):
        du = np.zeros(2); du[j] = eps_u
        B[:, j] = (f(x_op, u_op + du) - f(x_op, u_op - du)) / (2 * eps_u)
    return A, B

# ---------------------------------------------------------
# 2) Dobór PID (Wewnętrzna pętla - Pionowanie)
# ---------------------------------------------------------
def pid_from_pole_placement_for_psi(A, B, desired_poles):
    Bc = 0.5 * (B[:, 0] + B[:, 1])
    Apsi = A[np.ix_([2, 3], [2, 3])]
    Bpsi = Bc[[2, 3]].reshape(2, 1)

    # Augmentacja o całkę (z poprawionym znakiem dla spójności)
    Aaug = np.array([
        [Apsi[0, 0], Apsi[0, 1], 0.0],
        [Apsi[1, 0], Apsi[1, 1], 0.0],
        [1.0,        0.0,        0.0], 
    ])
    Baug = np.vstack([Bpsi, [[0.0]]])
    
    placed = place_poles(Aaug, Baug, desired_poles)
    K = placed.gain_matrix.reshape(-1)
    return float(K[0]), float(K[1]), float(K[2])

# -----------------------------
# 3) Logika Kaskady (Nowa funkcja)
# -----------------------------
def cascade_step(t, x, x_ref, pid_inner, kp_pos, kd_pos):
    """
    Pętla zewnętrzna: steruje pozycją kół (theta -> x[0])
    Wylicza psi_ref dla pętli wewnętrznej.
    """
    theta_error = x_ref[0] - x[0]
    theta_dot = x[1]
    
    # Wyliczamy o ile musimy się pochylić, żeby wrócić na pozycję
    # Używamy abs(), żeby dopasować znak do Twojego kp_inner
    sign_fix = np.sign(pid_inner.g.kp)
    psi_ref = sign_fix * (kp_pos * theta_error - kd_pos * theta_dot)
    
    # Saturacja nachylenia - max 10 stopni (0.17 rad)
    psi_ref = np.clip(psi_ref, -0.17, 0.17)
    
    # Pętla wewnętrzna (Twoja klasa)
    # Tworzymy wektor referencyjny z nowym psi_ref
    x_ref_inner = np.copy(x_ref)
    x_ref_inner[2] = psi_ref 
    
    return pid_inner(t, x, x_ref_inner)

def cascade_step(t, x, x_ref, pid_inner, kp_pos, kd_pos):
    # 1. Obliczamy błąd pozycji (theta to x[0])
    theta_error = x_ref[0] - x[0]
    theta_dot = x[1]
    
    # 2. KLUCZOWA ZMIANA: Odwracamy znak. 
    # Jeśli robot uciekał w jedną stronę, musimy zmienić znak wzmocnienia.
    # W większości modeli Segwaya: psi_ref = -(KP*err + KD*vel)
    # zestaw z różnymi znakami wzmocnień xd (coś jest nie tak) chyba jest dobrze, bo błąd pcohodnej to: e_theta_dot == theta_dot_ref- theta_dot = - theta_dot
    psi_ref = (kp_pos * theta_error - kd_pos * theta_dot)
    # # eksperymenty
    # psi_ref = -(kp_pos * theta_error + kd_pos * theta_dot)
    
    # 3. Ograniczenie (Saturacja)
    psi_ref = np.clip(psi_ref, -0.17, 0.17)
    
    # 4. Przygotowanie wektora referencji dla PIDController
    # Twoja klasa PID szuka wartości zadanej na indeksie self.psi_idx (domyślnie 2)
    x_ref_inner = np.zeros(6)
    x_ref_inner[pid_inner.psi_idx] = psi_ref
    
    return pid_inner.compute(t, x, x_ref_inner)

# -----------------------------
# 4) Main / Symulacja
# -----------------------------
if __name__ == "__main__":
    model = UnTransNonlinearModel()
    A, B = linearize_fd(model, np.zeros(6), np.zeros(2))

    # Nastawy pętli wewnętrznej (Pionowanie)
    inner_poles = (-2.5, -2.6, -2.7)
    kp, kd, ki = pid_from_pole_placement_for_psi(A, B, inner_poles)
    
    print(f"Inner PID (Psi): Kp={kp:.2f}, Ki={ki:.2f}, Kd={kd:.2f}")

    pid_inner = PIDController(
        gains=PIDGains(kp=kp, ki=ki, kd=kd),
        dt=0.001, psi_idx=2, psi_dot_idx=3,
        u_min=-12.0, u_max=12.0
    )

    # Parametry pętli zewnętrznej (Pozycja) - dobierane eksperymentalnie
    KP_POS = 0.00002  
    KD_POS = 0.000001

    # zestaw z różnymi znakami wzmocnień xd (coś jest nie tak) jest dobrze, bo błąd pcohodnej to: e_theta_dot == theta_dot_ref- theta_dot = - theta_dot
    KP_POS = 0.002  
    KD_POS = 0.001 # bardzo dobrze


    # Warunek początkowy: robot lekko wychylony
    x0 = np.array([0, 0, 5*np.pi/180, 0, 0, 0])
    x_target = np.array([0, 0, 0, 0, 0, 0]) # Cel: theta=0, psi=0

    def system_dynamics(t, x):
        u = cascade_step(t, x, x_target, pid_inner, KP_POS, KD_POS)
        return model.derivatives(t, x, u[0], u[1])

    sol = solve_ivp(system_dynamics, (0, 150), x0, max_step=0.005)

    # Wykresy
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 8))
    ax1.plot(sol.t, sol.y[0], label="Kąt kół theta (pozycja)")
    ax1.set_ylabel("rad")
    ax1.legend(); ax1.grid(True)
    
    ax2.plot(sol.t, sol.y[2], color='orange', label="Kąt nachylenia psi (pion)")
    ax2.set_ylabel("rad")
    ax2.set_xlabel("Czas [s]")
    ax2.legend(); ax2.grid(True)
    
    plt.suptitle("Sterowanie Kaskadowe: Pozycja + Pionowanie")
    plt.show()