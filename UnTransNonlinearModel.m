classdef UnTransNonlinearModel < handle
    properties
        p % Parametry modelu
    end
    
    methods
        function obj = UnTransNonlinearModel(params)
            if nargin < 1 || isempty(params)
                obj.p = obj.default_params();
            else
                obj.p = params;
            end
        end
        
        function p = default_params(obj)
            % Masa i geometria
            p.m = 0.32;             % masa jednego koła [kg]
            p.M = 5.41;             % masa "nadwozia" [kg]
            p.R = 0.15 / 2.0;       % promień koła [m]
            p.W = 0.40;             % rozstaw kół [m] 
            p.Jphi = 0.0484;        % bezwładność osi obrotu w planie [kg*m^2]
            p.Jm = 0.00119;         % bezwładność silnik+przekładnia [kg*m^2]

            % Silnik i tarcie
            p.Rm = 1.0;             % rezystancja uzwojenia [ohm]
            p.Kt = 0.025;           % stała momentu [Nm/A]
            p.Ke = 0.025;           % stała SEM [Vs/rad]
            p.fm = 0.00024;         % współczynnik tarcia

            % Grawitacja
            p.g = 9.81;             % [m/s^2]
        end
        
        function [d1, d2] = get_d1_d2(obj)
            d1 = obj.p.Kt / obj.p.Rm;
            d2 = (obj.p.Kt * obj.p.Ke) / obj.p.Rm + obj.p.fm;
        end
        
        function x_dot = derivatives(obj, t, x, u_l, u_r)
            % x = [theta, theta_dot, psi, psi_dot, phi, phi_dot]
            theta_dot = x(2);
            psi       = x(3);
            psi_dot   = x(4);
            phi_dot   = x(6);

            p = obj.p;
            [d1, d2] = obj.get_d1_d2();

            % Współczynniki (4.31–4.32)
            c1 = ((2.0 * p.m + p.M) * (p.R^2) + 2.0 * p.Jw + 2.0 * p.Jm);
            c2 = (p.M * (p.L^2) + p.Jpsi + 2.0 * p.Jm);
            f1 = (p.M * p.R * p.L * cos(psi) - 2.0 * p.Jm);
            f2 = (p.M * p.R * p.L * (psi_dot^2) * sin(psi));
            f3 = (p.M * (p.L^2) * (phi_dot^2) * sin(psi) * cos(psi) + p.M * p.g * p.L * sin(psi));

            denom = (c1 * c2 - f1 * f1);
            if abs(denom) < 1e-9
                denom = sign(denom) * 1e-9;
                if denom == 0, denom = 1e-9; end
            end

            % Siły uogólnione (4.36–4.38)
            F_theta = d1 * (u_r + u_l) + 2.0 * d2 * (psi_dot - theta_dot);
            F_psi   = -d1 * (u_r + u_l) + 2.0 * d2 * (theta_dot - psi_dot);
            F_phi   = -d1 * (p.W / (2.0 * p.R)) * (u_r - u_l) - ((p.W^2) / (2.0 * (p.R^2))) * d2 * phi_dot;

            % Przyspieszenia (4.33–4.35)
            theta_ddot = (c2 * (F_theta + f2) - f1 * (F_psi + f3)) / denom;
            psi_ddot   = (c1 * (F_psi + f3) - f1 * (F_theta + f2)) / denom;

            denom_phi = (0.5 * p.m * (p.W^2) + ...
                         p.M * (p.L^2) * (sin(psi)^2) + ...
                         (p.W^2) / (2.0 * (p.R^2)) * (p.Jw + p.Jm) + ...
                         p.Jphi);

            if abs(denom_phi) < 1e-12
                denom_phi = 1e-12;
            end

            phi_ddot = (F_phi - 2.0 * p.M * (p.L^2) * phi_dot * psi_dot * sin(psi) * cos(psi)) / denom_phi;

            % Wynik jako wektor kolumnowy
            x_dot = [theta_dot; theta_ddot; psi_dot; psi_ddot; phi_dot; phi_ddot];
        end
        
        function x_next = step_euler(obj, x, u_l, u_r, dt, t)
            if nargin < 6, t = 0; end
            dx = obj.derivatives(t, x, u_l, u_r);
            x_next = x + dt * dx;
        end
    end
end