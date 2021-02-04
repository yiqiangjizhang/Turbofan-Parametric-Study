%% Optimization Solver
% Second method for finding optimal parameters
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
    Tt4 = 1450;         % Inlet turbine temperature [K]
    
%% 2. Solver with design parameters
    alpha_d = 8;
    pi_f_d = 1.5;
    pi_lpc_d = 4;
    pi_hpc_d = 8.5;
    [f,T,P,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp] = ...
            solver(M0,altitude,Tt4,alpha_d,pi_f_d,pi_lpc_d,pi_hpc_d);
    F_spec_d=F_spec_total;
    I_sp_d=I_sp;
    
%%  2. Solver with alpha as a variable
    
    alpha_vec = 4:0.1:14.5;
    pi_f = 1.5;
    pi_lpc = 4;
    pi_hpc = 8.5;
    F_spec_vector = zeros(1,length( alpha_vec));
    I_sp_vector = zeros(1,length(alpha_vec));
    for i = 1:length(alpha_vec)

        [f,T,P,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp] = ...
            solver(M0,altitude,Tt4,alpha_vec(i),pi_f,pi_lpc,pi_hpc);

         F_spec_vector(i) = F_spec_total;
         I_sp_vector(i) = I_sp;
    end
    max1=max(F_spec_vector);
    position1=find(F_spec_vector==max1);
    max2=max(I_sp_vector);
    position2=find(I_sp_vector==max2);
    alpha_opt=0.1*alpha_vec(position1)+0.9*alpha_vec(position2);
    
%%  3. Solver with pi_f as a variable

    pi_f_vec = 1.25:0.01:2;
    pi_lpc = 4;
    pi_hpc = 8.5;
    F_spec_vector = zeros(1,length(pi_f_vec));
    I_sp_vector = zeros(1,length(pi_f_vec));
    for i = 1:length(pi_f_vec)
        [f,T,P,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp] = ...
            solver(M0,altitude,Tt4,alpha_opt,pi_f_vec(i),pi_lpc,pi_hpc);

         F_spec_vector(i) = F_spec_total;
         I_sp_vector(i) = I_sp;

         if imag(F_spec_vector(i))~=0
             F_spec_vector(i)=0;
         end
         if imag(I_sp_vector(i))~=0
             I_sp_vector(i)=0;
         end    
    end
    max1=max(F_spec_vector);
    position1=find(F_spec_vector==max1);
    max2=max(I_sp_vector);
    position2=find(I_sp_vector==max2);
    pi_f_opt=0.1*pi_f_vec(position1)+0.9*pi_f_vec(position2); %si creemos que el thrust es mas importante
                                                             %podemos hacer una media ponderada
%%  4. Solver with pi_hpc as a variable

    pi_hpc_vec = 1.3:0.01:12;
    pi_lpc = 4;
    F_spec_vector = zeros(1,length( pi_hpc_vec));
    I_sp_vector = zeros(1,length(pi_hpc_vec));
    for i = 1:length(pi_hpc_vec)

        [f,T,P,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp] = ...
            solver(M0,altitude,Tt4,alpha_opt,pi_f_opt,pi_lpc,pi_hpc_vec(i));

         F_spec_vector(i) = F_spec_total;
         I_sp_vector(i) = I_sp;

    end
    max1=max(F_spec_vector);
    position1=find(F_spec_vector==max1);
    max2=max(I_sp_vector);
    position2=find(I_sp_vector==max2);
    pi_hpc_opt=0.1*pi_hpc_vec(position1)+0.9*pi_hpc_vec(position2);
% figure(1)
% plot(pi_hpc_vec,F_spec_vector,'r')
% title('Fuerza espec�fica vs pi_hpc')
% 
% figure(2)
% plot(pi_hpc_vec,I_sp_vector,'b')
% title('Impulso espec�fico vs pi_hpc')

%%  4. Solver with pi_lpc as a variable

    pi_lpc_vec = 1.3:0.01:6;
    F_spec_vector = zeros(1,length( pi_lpc_vec));
    I_sp_vector = zeros(1,length(pi_lpc_vec));
    for i = 1:length(pi_lpc_vec)

        [f,T,P,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp] = ...
            solver(M0,altitude,Tt4,alpha_opt,pi_f_opt,pi_lpc_vec(i),pi_hpc_opt);

         F_spec_vector(i) = F_spec_total;
         I_sp_vector(i) = I_sp;

    end
    max1=max(F_spec_vector);
    F_spec_opt=max1;
    position1=find(F_spec_vector==max1);
    max2=max(I_sp_vector);
    I_sp_opt=max2;
    position2=find(I_sp_vector==max2);
    pi_lpc_opt=0.1*pi_lpc_vec(position1)+0.9*pi_lpc_vec(position2);
    
 %% GAIN OF ADIMENSIONAL THRUST AND SPECIFIC IMPULSE
    [f,T,P,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp] = ...
            solver(M0,altitude,Tt4,alpha_opt,pi_f_opt,pi_lpc_opt,pi_hpc_opt);
    F_spec_opt=F_spec_total;
    I_sp_opt=I_sp;
    % Comparision
    Case = {'Design';'Optimized'};
    F_spec_vecopt=[F_spec_d;F_spec_opt];
    I_spec_vecopt=[I_sp_d;I_sp_opt];
    alpha_vecopt=[alpha_d;alpha_opt];
    pi_f_vecopt=[pi_f_d;pi_f_opt];
    pi_hpc_vecopt=[pi_hpc_d;pi_hpc_opt];
    pi_lpc_vecopt=[pi_lpc_d;pi_lpc_opt];
    T=table(Case,F_spec_vecopt, I_spec_vecopt,alpha_vecopt,pi_f_vecopt,pi_hpc_vecopt,pi_lpc_vecopt)

    
    fprintf("Optimal values\n");
    fprintf("%15s = %.4f\n", "alpha_opt", alpha_opt);
    fprintf("%15s = %.4f\n", "pi_f_opt", pi_f_opt);
    fprintf("%15s = %.4f\n", "pi_hpc_opt", pi_hpc_opt);
    fprintf("%15s = %.4f\n", "pi_lpc_opt", pi_lpc_opt);
    
%  %% 5. Propulsive efficiency and exhaust velocities
%     %Valors optims
%     M0 = 0.85;
%     altitude = 11000;
%     Tt4 = 1450;
%     alpha_opt=12;
%     pi_lpc_opt=2.8;
%     pi_hpc_opt=15;
%     pi_f_vec = 1:0.01:2.2;
%     prop_eff = zeros(1,length(pi_f_vec));
%     U9_vec = zeros(1,length(pi_f_vec));
%     U1_9_vec = zeros(1,length(pi_f_vec));
%     for i = 1:length(pi_f_vec)
% 
%         [f,T,P,U0,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp,mfp_9,mfp1_9,mfp_8,mfp1_8] = ...
%             solver_matched(M0,altitude,Tt4,alpha_opt,pi_f_vec(i),pi_lpc_opt,pi_hpc_opt);
%          if isreal(U9)
%              U9_vec(i) = U9;
%          end
%          if isreal(U1_9)
%              U1_9_vec(i)=U1_9;
%          end
%          if isreal(U9)&& isreal(U1_9)
%              prop_eff(i) = 2*F_spec_total*U0/((1+f)*U9^2+alpha_opt*U1_9^2-(1+alpha_opt)*U0^2);
%          end
%          
%     end
% 
% %Plots  pi_f_vec = 1:0.01:2.2;
%     prop_eff = zeros(1,length(pi_f_vec));
%     U9_vec = zeros(1,length(pi_f_vec));
%     U1_9_vec = zeros(1,length(pi_f_vec));
%     for i = 1:length(pi_f_vec)
% 
%         [f,T,P,U0,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp,mfp_9,mfp1_9,mfp_8,mfp1_8] = ...
%             solver_matched(M0,altitude,Tt4,alpha_opt,pi_f_vec(i),pi_lpc_opt,pi_hpc_opt);
%          if isreal(U9)
%              U9_vec(i) = U9;
%          end
%          if isreal(U1_9)
%              U1_9_vec(i)=U1_9;
%          end
%          if isreal(U9)&& isreal(U1_9)
%              prop_eff(i) = 2*F_spec_total*U0/((1+f)*U9^2+alpha_opt*U1_9^2-(1+alpha_opt)*U0^2);
%          end
%          
%     end
% 
% figure(1)
% hold on
% title('\textbf{Propulsive efficiency vs. Fan Pressure Ratio}');
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% plot(pi_f_vec, prop_eff, 'r', 'LineWidth', 1);
% xlim([1 1.7]);
% xlabel("Fan Ratio");
% ylabel("Propulsive efficiency");
% set(gcf,'units','centimeters','position',[1,1,18,15]);
% grid on;
% grid minor;
% box on;
% hold off;
% 
% figure(2)
% hold on
% title('\textbf{Exhaust velocities vs. Fan Pressure Ratio}');
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% plot(pi_f_vec, U9_vec, 'r', 'LineWidth', 1);
% plot(pi_f_vec, U1_9_vec, 'b', 'LineWidth', 1);
% 
% 
% xlim([1 2.2]);
% xlabel("Fan Pressure Ratio");
% ylabel("Exhaust velocity $\left( \mathrm{m} \cdot \mathrm{s} \right)$");
% set(gcf,'units','centimeters','position',[19,1,18,15]);
% grid on;
% grid minor;
% box on;
% legend('Primary flow exhaust velocity $\left( U_9 \right)$', 'Secondary flow exhaust velocity $\left( U_{19} \right)$', 'Location', 'northeast');
% hold off