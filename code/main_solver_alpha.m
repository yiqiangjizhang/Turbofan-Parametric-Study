%% main_solver_alpha
% Script to perform sensitivity analysis of specific thrust and specific 
% impulse with alpha
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

alpha = 3.5:0.1:14.5; % Bypass ratio
pi_f = 1.5;         % Fan pressure ratio
pi_lpc = 4;         % LPC pressure ratio
pi_hpc = 8.5;       % HPC pressure ratio
Tt4 = 1450;         % Inlet turbine temperature [K]

F_spec_vector = zeros(1,length(alpha)); % Specific thrust vector
I_sp_vector = zeros(1,length(alpha));   % Specific impulse vector

%%  2. Solver
% 2.1. Solve for a certain range of alpha
for i = 1:length(alpha)
    % Solve general turbofan
    [f,T,P,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp] = ...
        solver(M0,altitude,Tt4,alpha(i),pi_f,pi_lpc,pi_hpc);
    % Obtain parameters
    F_spec_vector(i) = F_spec_total;
    I_sp_vector(i) = I_sp; 
end
 
%% 3. PLOTS
save = false;

% 3.1. Specific thrust vs alpha
figure(1);
hold on;
title('\textbf{Specific Thrust vs. Bypass Ratio}');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(alpha, F_spec_vector, 'r', 'LineWidth', 1);
ylim([950 1450]);
xlabel("Bypass Ratio");
ylabel("Specific Thrust $\left( \mathrm{m} \cdot \mathrm{s}^{-1} \right)$");
set(gcf,'units','centimeters','position',[1,1,18,15]);
grid on;
grid minor;
box on;
if save == true
    saveas(gcf, 'plots/plot_spec_thrust_alpha.svg')
end
hold off;

% 3.2. Specific impulse vs alpha
figure(2);
hold on;
title('\textbf{Specific Impulse vs. Bypass Ratio}');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(alpha, I_sp_vector, 'b', 'LineWidth', 1);
ylim([3800 5800]);
xlabel("Bypass Ratio");
ylabel("Specific Impulse $\left( \mathrm{s} \right)$");
set(gcf,'units','centimeters','position',[19,1,18,15]);
grid on;
grid minor;
box on;
if save == true
    saveas(gcf, 'plots/plot_spec_impulse_alpha.svg')
end
hold off;





