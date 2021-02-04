%% main
% Computes thermal, propulsive and global turbofan efficiencies

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

%% 5. Propulsive Efficiency

pi_f_vec = 1:0.01:2;
prop_eff = zeros(1,length(pi_f_vec));
U9_vec = zeros(1,length(pi_f_vec));
U1_9_vec = zeros(1,length(pi_f_vec));

for i = 1:length(pi_f_vec)

    [f,T,P,U0,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp,mfp_9,mfp1_9,mfp_8,mfp1_8] = ...
        solver_matched(M0,altitude,Tt4,alpha,pi_f_vec(i),pi_lpc,pi_hpc);

    if pi_f_vec(i) == 1.5
        U0_e = U0;
        f_e = f;
        F_spec_total_e = F_spec_total;
        U_9_e = U9;
        U1_9_e = U1_9;
    end

    if isreal(U9)
        U9_vec(i) = U9;
    end

    if isreal(U1_9)
        U1_9_vec(i)=U1_9;
    end

    if isreal(U9)&& isreal(U1_9)
        prop_eff(i) = 2*F_spec_total*U0/((1+f)*U9^2+alpha*U1_9^2-(1+alpha)*U0^2);
    end

end

figure(1);
plot(pi_f_vec,prop_eff,'r');
title('Propulsive efficiency vs pi_f');

figure(2);
plot(pi_f_vec,U9_vec,'b',pi_f_vec,U1_9_vec,'g');
title('Exhaust velocities vs pi_f');
legend('U_9','U_{19}');

thermal_eff_1 = ((1+f_e)*U_9_e^2+alpha*U1_9_e^2-(1+alpha)*U0_e^2)/(2*f_e*h);
prop_eff_1 = 2*F_spec_total_e*U0_e/((1+f_e)*U_9_e^2+alpha*U1_9_e^2-(1+alpha)*U0_e^2);
glob_ef_1 = thermal_eff_1*prop_eff_1;

fprintf("Efficiencies\n");
fprintf("%15s = %.4f\n", "thermal_eff_1", thermal_eff_1);
fprintf("%15s = %.4f\n", "prop_eff_1", prop_eff_1);
fprintf("%15s = %.4f\n", "glob_ef_1", glob_ef_1);
