import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parametry z robota ETH (z Table 2, page 10 w mpc_nmpc_pid_lqr_genialne.pdf)
d = 0.1      # odległość między kołami [m]
l = 0.0275   # odległość osi kół do CoM [m]
r = 0.04     # promień koła [m]
m_B = 0.43   # masa ciała [kg]
m_W = 0.02   # masa koła [kg]
J = 1.92e-4  # MOI koła wzdłuż osi [kg m²]
g = 9.81     # grawitacja [m/s²]
K = 1e-4     # MOI koła w pionie [kg m²]
I2 = 1.49e-3 # MOI ciała wzdłuż B [kg m²]

# Uproszczone stałe (dla 2D)
m_o = m_B + 2 * m_W + 2 * J / r**2
I_o = I2 + m_B * l**2
a = m_B * l

# Nieliniowa dynamika 2D (stan: [x, dx, theta, dtheta], input T = T_L + T_R)
def dynamics(state, t, T):
    x, dx, theta, dtheta = state
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    denom = m_o * I_o - a**2 * cos_theta**2
    if denom == 0: denom = 1e-6  # unikaj dzielenia przez 0
    
    ddx = (a * I_o * dtheta**2 * sin_theta - a**2 * g * sin_theta * cos_theta + T * (I_o / r + a * cos_theta)) / denom
    ddtheta = (-a**2 * dtheta**2 * sin_theta * cos_theta + a * m_o * g * sin_theta - T * (m_o + a / r * cos_theta)) / denom
    
    return np.array([dx, ddx, dtheta, ddtheta])

# SMC (Sliding Mode Control) - z artykułu segway_genetic_algorithm.pdf
def smc(state, c1, c2, epsilon, theta_d=0, dtheta_d=0, dx_d=0):
    dx, theta, dtheta = state[1], state[2], state[3]
    e = np.array([dx - dx_d, theta - theta_d, dtheta - dtheta_d])
    C = np.array([c1, c2, 1.0])  # c3=1 jak w artykule
    s = np.dot(C, e)
    T = -epsilon * np.sign(s) * (m_o + a / r)  # uproszczone prawo sterujące (torque, dostosowane do stabilizacji)
    return T

# PID (kaskadowy: zewnętrzny dla x, wewnętrzny dla theta)
def pid_cascaded(state, integral_theta, integral_x, kp_theta, ki_theta, kd_theta, kp_x, ki_x, kd_x, dt, x_d=0, theta_d=0):
    x, dx, theta, dtheta = state
    # Wewnętrzny PID dla theta (referencja z zewnętrznego)
    theta_ref = kp_x * (x_d - x) + ki_x * integral_x + kd_x * (-dx)  # ref theta z błędu x
    e_theta = theta_ref - theta
    integral_theta += e_theta * dt
    de_theta = -dtheta
    T = kp_theta * e_theta + ki_theta * integral_theta + kd_theta * de_theta # z minusem lepiej, ale testuj oba
    
    integral_x += (x_d - x) * dt
    return T, integral_theta, integral_x

# Switching Control
def switching_controller(state, alpha, c1, c2, epsilon, kp_theta, ki_theta, kd_theta, kp_x, ki_x, kd_x, integral_theta, integral_x, dt, x_d=0, theta_d=0):
    theta = state[2]
    if abs(theta) > alpha:
        T = smc(state, c1, c2, epsilon, theta_d, dtheta_d=0, dx_d=0)
    else:
        T, integral_theta, integral_x = pid_cascaded(state, integral_theta, integral_x, kp_theta, ki_theta, kd_theta, kp_x, ki_x, kd_x, dt, x_d, theta_d)
    return T, integral_theta, integral_x

# Symulacja odpowiedzi (dla fitness, zwraca fitness)
def simulate(params, initial_state=np.array([0, 0, np.pi/4, 0]), t_span=3.0, dt=0.004, test_type='angle_recovery'):
    c1, c2, epsilon, kp_theta, ki_theta, kd_theta, kp_x, ki_x, kd_x, alpha = params
    t = np.arange(0, t_span, dt)
    state = initial_state
    states = [state]
    integral_theta = 0
    integral_x = 0
    Ts = []
    for i in range(1, len(t)):
        T, integral_theta, integral_x = switching_controller(state, alpha, c1, c2, epsilon, kp_theta, ki_theta, kd_theta, kp_x, ki_x, kd_x, integral_theta, integral_x, dt)
        state = odeint(dynamics, state, [t[i-1], t[i]], args=(T,))[-1]
        states.append(state)
        Ts.append(T)
    states = np.array(states)
    Ts = np.array(Ts)
    
    # Fitness z artykułu (eq 7 dla angle recovery)
    theta = states[:, 2]
    theta_max = np.max(np.abs(theta))
    theta_tot = np.trapz(np.abs(theta), t)
    d_trav = np.max(np.abs(states[:, 0]))
    u_max = np.max(np.abs(states[:, 1]))
    E = np.trapz(Ts**2, t[:-1])  # proxy energii (T^2)
    tu = np.sum(np.abs(theta) < 0.05) * dt  # czas upright
    fitness = 50 * tu / (E + 1e-6) * (100 * theta_max + theta_tot + u_max + d_trav)  # max fitness
    return -fitness  # dla min w GA (wyższa fitness lepsza)

# Genetic Algorithm (prosty, z numpy)
def genetic_algorithm(pop_size=50, gens=20, test_type='angle_recovery'):
    # Granice parametrów (z artykułu: c1,c2 ±1, ε>0, PID>0, alpha 0-pi/6)
    bounds = [(-1,1), (-1,1), (0,10), (10,100), (0,10), (10,100), (0,50), (0,10), (10,100), (0, np.pi/6)]
    pop = np.random.uniform([b[0] for b in bounds], [b[1] for b in bounds], (pop_size, 10))
    best_fitness = []
    for gen in range(gens):
        fitness = np.array([simulate(p, test_type=test_type) for p in pop])
        best_fitness.append(-np.min(fitness))  # track best
        idx = np.argsort(fitness)[:pop_size//2]  # select top half (min negative = max fitness)
        parents = pop[idx]
        offspring = []
        for _ in range(pop_size - len(parents)):
            p1, p2 = parents[np.random.choice(len(parents), 2)]
            child = (p1 + p2) / 2 + np.random.normal(0, 0.05, 10)  # crossover + mutation
            child = np.clip(child, [b[0] for b in bounds], [b[1] for b in bounds])
            offspring.append(child)
        pop = np.vstack((parents, offspring))
    best_idx = np.argmin(fitness)
    return pop[best_idx], best_fitness

# Funkcja do symulacji z plotem (bez fitness)
def run_simulation(params, initial_state=np.array([0, 0, np.pi/4, 0]), t_span=3.0, dt=0.004):
    c1, c2, epsilon, kp_theta, ki_theta, kd_theta, kp_x, ki_x, kd_x, alpha = params
    t = np.arange(0, t_span, dt)
    state = initial_state
    states = [state]
    integral_theta = 0
    integral_x = 0
    Ts = []
    for i in range(1, len(t)):
        T, integral_theta, integral_x = switching_controller(state, alpha, c1, c2, epsilon, kp_theta, ki_theta, kd_theta, kp_x, ki_x, kd_x, integral_theta, integral_x, dt)
        state = odeint(dynamics, state, [t[i-1], t[i]], args=(T,))[-1]
        states.append(state)
        Ts.append(T)
    states = np.array(states)
    Ts = np.array(Ts)
    return states, t, Ts

# Przykład użycia
# 1. Uruchom GA
best_params, fitness_history = genetic_algorithm(pop_size=20, gens=10)  # mniejszy dla testu
print("Optymalne parametry:", best_params)
print("Historia fitness:", fitness_history)

# 2. Symuluj z najlepszymi
states, t, Ts = run_simulation(best_params, initial_state=np.array([0, 0, np.deg2rad(45), 0]), t_span=5.0, dt=0.01)

# 3. Plot
plt.figure(figsize=(12, 8))
plt.subplot(3,1,1)
plt.plot(t, np.rad2deg(states[:,2]), label='Kąt θ [deg]')
plt.axhline(0, color='r', linestyle='--')
plt.legend()
plt.title('Symulacja z GA-optymalizowanymi parametrami')
plt.subplot(3,1,2)
plt.plot(t, states[:,0], label='Pozycja x [m]')
plt.legend()
plt.subplot(3,1,3)
plt.plot(t[:-1], Ts, label='Torque T [Nm]')
plt.legend()
plt.tight_layout()
plt.show()