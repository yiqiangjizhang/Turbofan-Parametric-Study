%% main_solver_alpha
% Solves the turbojet engine given the input values
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
Tt4 = 1450;         % Inlet turbine temperature [K]

%%  2. Solver for design parameters
%Valors disseny: alpha = 8; pi_lpc = 4; pi_hpc = 8; pi_f = 1.5
%Valors optims: alpha = 12; pi_lpc = 2.8; pi_hpc = 15; pi_f = 1.5
alpha=0;
pi_f=0;
pi_lpc = 5;  
pi_hpc = 10;
[f,T,P,U9,U1_9,M9,M1_9,F_spec_total,c_s,I_sp] = ...
    solver_turbojet(M0,altitude,Tt4,alpha,pi_f,pi_lpc,pi_hpc);
F_spec_total
I_sp
%%  3. Solver for pi_lpc
%Valors disseny: alpha = 8; pi_lpc = 4; pi_hpc = 8; pi_f = 1.5
alpha=0;
pi_f=0;
pi_lpc = 1.3:0.01:6.25;  
pi_hpc = 8;
F_spec_vector = zeros(1,length(pi_lpc)); % Specific thrust vector
I_sp_vector = zeros(1,length(pi_lpc));   % Specific impulse vector
% 2.1. Solve for a certain range of pi_hpc
for i = 1:length(pi_lpc)
    % Solve general turbofan    
    [f,T,P,U9,U1_9,M9,M1_9,F_spec_total,c_s,I_sp] = ...
    solver_turbojet(M0,altitude,Tt4,alpha,pi_f,pi_lpc(i),pi_hpc);
    % Obtain parameters
    F_spec_vector(i) = F_spec_total;
    I_sp_vector(i) = I_sp; 
end
%% 4. PLOTS for pi_lpc
save = false;

% 3.1. Specific thrust vs LPC Pressure Ratio
figure(1);
hold on;
title('\textbf{Specific Thrust vs. LPC Pressure ratio}');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(pi_lpc, F_spec_vector, 'r', 'LineWidth', 1);
%ylim([1180 1320]);
xlabel("LPC Pressure Ratio");
ylabel("Specific Thrust $\left( \mathrm{m} \cdot \mathrm{s}^{-1} \right)$");
set(gcf,'units','points','position',[50,50,550,400]);
grid on;
grid minor;
box on;
if save == true
    saveas(gcf, 'plots/plot_spec_thrust_pi_lpc.svg')
end
hold off;

% 3.2. Specific impulse vs LPC Pressure Ratio
figure(2);
hold on;
title('\textbf{Specific Impulse vs. LPC Pressure Ratio}');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(pi_lpc, I_sp_vector, 'b', 'LineWidth', 1);
xlabel("LPC Pressure Ratio");
ylabel("Specific Impulse $\left( \mathrm{s} \right)$");
set(gcf,'units','points','position',[600,50,550,400]);
grid on;
grid minor;
box on;
if save == true
    saveas(gcf, 'plots/plot_spec_impulse_pi_lpc.svg')
end
hold off;

%All plots
fig = figure(5);
subplot(2,2,1)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(pi_lpc, F_spec_vector, 'r', 'LineWidth', 1);
title('\textbf{Specific Thrust vs. LPC Pressure ratio}');
%ylim([1180 1320]);
xlabel("LPC Pressure Ratio");
ylabel("Specific Thrust $\left( \mathrm{m} \cdot \mathrm{s}^{-1} \right)$");
set(gcf,'units','centimeters','position',[1,1,18,15]);
grid on;
grid minor;
box on;

subplot(2,2,2)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(pi_lpc, I_sp_vector, 'b', 'LineWidth', 1);
title('\textbf{Specific Impulse vs. LPC Pressure Ratio}');
xlabel("LPC Pressure Ratio");
ylabel("Specific Impulse $\left( \mathrm{s} \right)$");
set(gcf,'units','centimeters','position',[1,1,18,15]);
grid on;
grid minor;
box on;
%%  5. Solver for pi_hpc
%Valors disseny: alpha = 8; pi_lpc = 4; pi_hpc = 8; pi_f = 1.5
alpha=0;
pi_f=0;
pi_lpc = 4;  
pi_hpc = 2:0.01:12;
F_spec_vector = zeros(1,length(pi_hpc)); % Specific thrust vector
I_sp_vector = zeros(1,length(pi_hpc));   % Specific impulse vector
% 2.1. Solve for a certain range of pi_hpc
for i = 1:length(pi_hpc)
    % Solve general turbofan    
    [f,T,P,U9,U1_9,M9,M1_9,F_spec_total,c_s,I_sp] = ...
    solver_turbojet(M0,altitude,Tt4,alpha,pi_f,pi_lpc,pi_hpc(i));
    % Obtain parameters
    F_spec_vector(i) = F_spec_total;
    I_sp_vector(i) = I_sp; 
end

%% 6. PLOTS for pi_hpc
save = false;

% 5.1. Specific thrust vs HPC Pressure Ratio
figure(3);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(pi_hpc, F_spec_vector, 'r', 'LineWidth', 1);
title('\textbf{Specific Thrust vs. HPC Pressure Ratio}');
xlim([1 13]);
ylim([680 780]);
xlabel("HPC Pressure Ratio");
ylabel("Specific Thrust $\left( \mathrm{m} \cdot \mathrm{s}^{-1} \right)$");
set(gcf,'units','centimeters','position',[1,1,18,15]);
grid on;
grid minor;
box on;
if save == true
    saveas(gcf, 'plots/plot_spec_thrust_pi_hpc_turbojet.svg')
end
hold off;

% 5.1. Specific impulse vs HPC Pressure Ratio
figure(4);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(pi_hpc, I_sp_vector, 'b', 'LineWidth', 1);
title('\textbf{Specific Impulse vs. HPC Pressure Ratio}');
xlim([1 13]);
%ylim([3800 5600]);
xlabel("HPC Pressure Ratio");
ylabel("Specific Impulse $\left( \mathrm{s} \right)$");
set(gcf,'units','centimeters','position',[1,1,18,15]);
grid on;
grid minor;
box on;
if save == true
    saveas(gcf, 'plots/plot_spec_impulse_pi_hpc_turbojet.svg')
end
hold off;

figure(5);
subplot(2,2,3)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(pi_hpc, F_spec_vector, 'r', 'LineWidth', 1);
title('\textbf{Specific Thrust vs. HPC Pressure Ratio}');
xlim([1 13]);
ylim([680 780]);
xlabel("HPC Pressure Ratio");
ylabel("Specific Thrust $\left( \mathrm{m} \cdot \mathrm{s}^{-1} \right)$");
set(gcf,'units','centimeters','position',[1,1,18,15]);
grid on;
grid minor;
box on;

subplot(2,2,4)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(pi_hpc, I_sp_vector, 'b', 'LineWidth', 1);
title('\textbf{Specific Impulse vs. HPC Pressure Ratio}');
xlim([1 13]);
%ylim([3800 5600]);
xlabel("HPC Pressure Ratio");
ylabel("Specific Impulse $\left( \mathrm{s} \right)$");
set(gcf,'units','centimeters','position',[1,1,18,15]);
grid on;
grid minor;
box on;
