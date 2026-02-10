%% UnTrans: stabilizacja + pozycja (PID-like MIMO)
%% UnTrans: Sterowanie PID-like MIMO
% Autor: Kamil Szady
% Data: 03-02-2026
% Opis: Skrypt realizuje linearyzację modelu nieliniowego oraz 
%       projektowanie regulatora metodą lokacji biegunów.
clear all;  % Usuwa zmienne z pamięci Workspace
close all;  % Zamyka wszystkie otwarte okna wykresów (poprawione)
clc;        % Czyści okno Command Window

% Inicjalizacja modelu (wymaga pliku UnTransNonlinearModel.m)
model = UnTransNonlinearModel(); 

% --- Parametry Symulacji ---
t_end = 30.0;
dt = 0.001;
u_min = -12.0;
u_max = 12.0;
z_limit = 50.0;

% --- Punkt Pracy i Referencja ---
x_op = zeros(6, 1);
u_op = zeros(2, 1);

s_ref = 0.0; 
theta_ref = s_ref / 0.1; % Przykładowy promień R=0.1, jeśli model go posiada
x_ref = [theta_ref; 0; 0; 0; 0; 0];

% --- 1) Linearyzacja (Finite Differences) ---
[A, B] = linearize_fd(model, x_op, u_op);

% --- 2) Projektowanie regulatora (Pole Placement) ---
y_indices = [1]; % Całkujemy tylko theta (indeks 1 w MATLAB)
ny = length(y_indices);

% Macierz C wybierająca wyjścia do całki
C = zeros(ny, 6);
for i = 1:ny
    C(i, y_indices(i)) = 1;
end

% Układ augmentowany
% [x_dot; z_dot] = [A 0; -C 0] * [x; z] + [B; 0] * u
Aaug = [A,           zeros(6, ny);
       -C,           zeros(ny, ny)];
Baug = [B;           zeros(ny, 2)];

% Bieguny (6 dla stanu + ny dla integratorów)
poles = [-10, -11, -12, -3, -3.5, -4, -0.5]; 

% Wyznaczanie wzmocnienia Kaug (u = -Kaug * [x; z])
% W MATLAB: place(A, B, p)
Kaug = place(Aaug, Baug, poles);

Kx = Kaug(:, 1:6);
Ki = Kaug(:, 7:end);

% --- 3) Symulacja Pętli Zamkniętej (RK4) ---
x0 = [0; 0; 5.5 * pi/180; 0; 0; 0]; % Warunek początkowy
z = zeros(ny, 1); % Stan integratorów

N = floor(t_end / dt) + 1;
t_arr = (0:N-1) * dt;
y_log = zeros(6, N);
u_log = zeros(2, N);

x = x0;
for k = 1:N
    y_log(:, k) = x;
    
    % Błąd i sterowanie
    e = x - x_ref;
    y_err = C * (x_ref - x);
    
    u_unsat = u_op - Kx * e - Ki * z;
    u_sat = max(u_min, min(u_max, u_unsat));
    u_log(:, k) = u_sat;
    
    % Anti-windup (Conditional Integration)
    if any(u_sat ~= u_unsat)
        dz = y_err * dt;
        for i = 1:ny
            du_i = -Ki(:, i) * dz(i);
            % Jeśli sterowanie jest nasycone i integracja by to pogłębiła - mrozimy
            worse_upper = (u_sat >= u_max - 1e-10) & (du_i > 0);
            worse_lower = (u_sat <= u_min + 1e-10) & (du_i < 0);
            if any(worse_upper | worse_lower)
                dz(i) = 0;
            end
        end
        z = z + dz;
    else
        z = z + y_err * dt;
    end
    z = max(-z_limit, min(z_limit, z));
    
    % Krok RK4
    if k < N
        x = rk4_step(model, t_arr(k), x, u_sat, dt);
    end
end

% --- 4) Wykresy ---
labels = {'\theta [rad]', 'd\theta/dt [rad/s]', '\psi [rad]', ...
          'd\psi/dt [rad/s]', '\phi [rad]', 'd\phi/dt [rad/s]'};
figure('Name', 'UnTrans - PID-like MIMO', 'NumberTitle', 'off');
for i = 1:6
    subplot(6, 1, i);
    plot(t_arr, y_log(i, :), 'LineWidth', 1.5);
    ylabel(labels{i});
    grid on;
end
xlabel('t [s]');

% --- Funkcje pomocnicze ---

function [A, B] = linearize_fd(model, x_op, u_op)
    eps_x = 1e-6;
    eps_u = 1e-6;
    A = zeros(6, 6);
    B = zeros(6, 2);
    
    %f = @(x, u) cell2mat(model.derivatives(0, x, u(1), u(2)));
    f = @(x, u) model.derivatives(0, x, u(1), u(2));
    
    for i = 1:6
        dx = zeros(6, 1); dx(i) = eps_x;
        A(:, i) = (f(x_op + dx, u_op) - f(x_op - dx, u_op)) / (2 * eps_x);
    end
    
    for j = 1:2
        du = zeros(2, 1); du(j) = eps_u;
        B(:, j) = (f(x_op, u_op + du) - f(x_op, u_op - du)) / (2 * eps_u);
    end
end

function x_next = rk4_step(model, t, x, u, dt)
    f = @(tt, xx) model.derivatives(tt, xx, u(1), u(2));
    
    k1 = f(t, x);
    k2 = f(t + dt/2, x + dt*k1/2);
    k3 = f(t + dt/2, x + dt*k2/2);
    k4 = f(t + dt, x + dt*k3);
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end