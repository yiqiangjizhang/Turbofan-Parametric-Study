%% Optimization Solver
% Method for finding optimal parameters

clc;
clear all;
close all;

%% GLOBAL PARAMETERS

global g ...
    pi_d eta_f eta_lpc eta_hpc pi_b eta_b eta_hpt eta_lpt eta_mh eta_ml pi_np pi_ns ...
    gamma_c gamma_t Rg h;

% Gravity Acceleration
g = 9.81;           % Earth gravity acceleration [m/s^2]

% Component eficiencies
pi_d = 0.98;        % Drag 
eta_f = 0.89;       % Fan isentropic efficiency
eta_lpc = 0.88;     % LPC isentropic efficiency [adim]
eta_hpc = 0.86;     % HPC isentropic efficiency [adim]
pi_b = 0.96;        % Pt4/Pt3 pressure ratio [adim]
eta_b = 0.99;       % Combustion efficiency [adim]
eta_hpt = 0.91;     % HPT isentropic efficiency [adim]
eta_lpt = 0.92;     % LPT isentropic efficiency [adim]
eta_mh = 0.993;     % Mechanical efficiency high section [adim]
eta_ml = 0.997;     % Mechanical efficiency low section [adim]
pi_np = 0.99;       % Primary nozzle efficiency [adim]
pi_ns = 0.99;       % Secondary nozzle efficiency [adim]

% Air and Gas properties
gamma_c = 1.4;      % Gamma compressor [adim]
gamma_t = 1.3;      % Gamma turbine [adim]
Rg = 287;           % Air constant [J/kgK]
h = 43e6;           % Combustion enthalpy [J/kg]

%% 1. Input data
% 1.1. Operational conditions
M0 = 0.85;          % Mach number [adim]
altitude = 11000;   % Altitude [m]

% 1.2. Design Parameters
alpha = 9:0.1:12;        % Bypass ratio
pi_f = 1.2:0.1:3;         % Fan pressure ratio
pi_lpc = 1.3:0.1:6;         % LPC pressure ratio
pi_hpc = 15:0.1:25;          % HPC pressure ratio
Tt4 = 1450;                 % Inlet turbine temperature [K]


%%  2. OPTIMIZER
x = 1;                    % Weight of F_spec_total
objective = intmin;         % Objective function
F_spec_total_opt = intmin;  % Optimal Specific thrust
I_sp_opt = intmin;          % Optimal Specific impulse
alpha_opt = intmin;         % Optimal bypass ratio
pi_f_opt = intmin;          % Optimal fan pressure ratio
pi_lpc_opt = intmin;        % Optimal LPC pressure ratio
pi_hpc_opt = intmin;        % Optimal HPC pressure ratio

fprintf("\nExpected time: %.2f %s\n", length(alpha)*length(pi_f)*length(pi_lpc)*length(pi_hpc)*(6e-5)/60, "min");


tic
% Search for alpha
for i = 1:length(alpha)
    % Search for pi_f
    for j = 1:length(pi_f)
        % Search for pi_lpc
        for k = 1:length(pi_lpc)
            % Search for pi_hpc
            for t = 1:length(pi_hpc)
                % Solve general turbofan
                [f,T,P,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp] = ...
                    solver(M0,altitude,Tt4,alpha(i),pi_f(j),pi_lpc(k),pi_hpc(t));
                % Check if output arguments are all real numbers
                if isreal(f) && isreal(T) && isreal(P) && isreal(U9) && isreal(U1_9) ...
                        && isreal(M9) && isreal(M1_9) && isreal(F_spec_p) && isreal(F_spec_s) ...
                        && isreal(F_spec_total) && isreal(c_s) && isreal(I_sp)
                    % Compute new value of objective function
                    %objective_case = x*F_spec_total + (1-x)*I_sp;
                    objective_case = F_spec_total*I_sp;
                    % Check if the new value is greater than the prior
                    if objective_case > objective
                        % Assign new values
                        objective = objective_case;
                        F_spec_total_opt = F_spec_total;
                        I_sp_opt = I_sp;
                        alpha_opt = alpha(i);
                        pi_f_opt = pi_f(j);
                        pi_lpc_opt = pi_lpc(k);
                        pi_hpc_opt = pi_hpc(t);                        
                    end
                end
            end
        end
    end
end
time = toc;

fprintf("\nOptimal solution after %d iterations\n", length(alpha)*length(pi_f)*length(pi_lpc)*length(pi_hpc));
fprintf("%20s = %.3f\n", "objective", objective);
fprintf("%20s = %.2f\n", "x", x);
fprintf("%20s = %.3f %s\n", "F_spec_total_opt", F_spec_total_opt, "m/s");
fprintf("%20s = %.3f %s\n", "I_sp_opt", I_sp_opt, "s");
fprintf("%20s = %.3f\n", "alpha_opt", alpha_opt);
fprintf("%20s = %.3f\n", "pi_f_opt", pi_f_opt);
fprintf("%20s = %.3f\n", "pi_lpc_opt", pi_lpc_opt);
fprintf("%20s = %.3f\n", "pi_hpc_opt", pi_hpc_opt);
fprintf("%20s = %.3f\n", "pi_global", pi_lpc_opt*pi_hpc_opt);
fprintf("%20s = %.3f %s\n", "time", time, "s");

fprintf("\n\n");

fprintf("%.2f\n", x);
fprintf("%.2f\n", F_spec_total_opt);
fprintf("%.2f\n", I_sp_opt);
fprintf("%.2f\n", alpha_opt);
fprintf("%.2f\n", pi_f_opt);
fprintf("%.2f\n", pi_lpc_opt);
fprintf("%.2f\n", pi_hpc_opt);
fprintf("%.2f\n", pi_lpc_opt*pi_hpc_opt);



%{
%% PLOTS
save = false;

% 3.1. Specific thrust vs Fan pressure ratio
figure(1);
hold on;
title('\textbf{Specific Thrust vs. Fan Pressure Ratio}');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(pi_f, F_spec_vector, 'r', 'LineWidth', 1);
xlim([1.1 2.1]);
xlabel("Fan Pressure Ratio");
ylabel("Specific Thrust $\left( \mathrm{m} \cdot \mathrm{s}^{-1} \right)$");
set(gcf,'units','points','position',[50,50,550,400]);
grid on;
grid minor;
box on;
if save == true
    saveas(gcf, 'plots/plot_spec_thrust_pi_f.svg')
end
hold off;

% 3.1. Specific impulse vs Fan pressure ratio
figure(2);
hold on;
title('\textbf{Specific Impulse vs. Fan Pressure Ratio}');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(pi_f, I_sp_vector, 'b', 'LineWidth', 1);
xlim([1.1 2.1]);
xlabel("Fan Pressure Ratio");
ylabel("Specific Thrust $\left( \mathrm{m} \cdot \mathrm{s}^{-1} \right)$");
set(gcf,'units','points','position',[600,50,550,400]);
grid on;
grid minor;
box on;
if save == true
    saveas(gcf, 'plots/plot_spec_impulse_pi_f.svg')
end
hold off;
%}
