%% main
% Solves the turbofan engine given the input values
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
alpha = 12;
pi_f = 1.5;
pi_lpc = 2.8;
pi_hpc= 15;
Tt4 = 1450;         % Inlet turbine temperature [K]

%%  2. Solver

[f,T,P,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp,I_sp_p,I_sp_s] = ...
    solver(M0,altitude,Tt4,alpha,pi_f,pi_lpc,pi_hpc);

%% 3. Solver with both nozzles matched

[f_m,T_m,P_m,U0,U9_m,U1_9_m,M9_m,M1_9_m,F_spec_p_m,F_spec_s_m,F_spec_total_m,c_s_m,I_sp_m,mfp_9,mfp1_9,mfp_8,mfp1_8] = ... 
     solver_matched(M0,altitude,Tt4,alpha,pi_f,pi_lpc,pi_hpc);

 
%% 4. Specific Impulse Gain

I_ps_Gain = 100*((I_sp_m-I_sp)/I_sp); % Specify impulse gain obtained with both nozzles matched
As_ratio = mfp1_8/mfp1_9;

fprintf("%15s = %.4f %s\n", "I_sp_m", I_sp_m, "s");
fprintf("%15s = %.4f %s\n", "I_sp", I_sp, "s");
fprintf("%15s = %.4f\n", "I_ps_Gain", I_ps_Gain);
fprintf("%15s = %.4f\n", "As_ratio", As_ratio);

