%%

function [f,T,P,U9,U1_9,M9,M1_9,F_spec_p,F_spec_s,F_spec_total,c_s,I_sp,I_sp_p,I_sp_s] = solver(M0,altitude,Tt4,alpha,pi_f,pi_lpc,pi_hpc)


global g ...
    pi_d eta_f eta_lpc eta_hpc pi_b eta_b eta_hpt eta_lpt eta_mh eta_ml pi_np pi_ns ...
    gamma_c gamma_t Rg h;

T = zeros(1,12);
P = zeros(1,10);

%% 1. Previous calculations

Cpc = Rg*gamma_c/(gamma_c-1);           % Specific Heat "compressor"
Cpt = Rg*gamma_t/(gamma_t-1);           % Specific Heat "turbine"

    % 1.1. Initial air conditions and paremeters
    T0 = 288.15-0.0065*altitude;            % Air temperature [K] %
    P0 = 101325*(T0/288.15)^5.256;          % Air pressure [Pa]
    theta_0 = 1+M0^2*(gamma_c-1)/2;         % Stagnation to static temperature ratio [adim]
    delta_0 = theta_0^(gamma_c/(gamma_c-1));% 
    Tt0 = T0*theta_0;                       %
    Pt0 = P0*delta_0;                       %
    U0 = M0*sqrt(gamma_c*Rg*T0);            % External air velocity [m/s]
    
%% 2. Solver

    % 3.1. Inlet
    Tt2 = Tt0;          % Supposed isoentropyc proces while inlet flow
    Pt2 = pi_d*Pt0;     % 

    % 3.2. SECONDARY FLOW
    
    % 3.2.1. Fan's Compression 
    tau_f = 1+((pi_f^((gamma_c-1)/gamma_c)-1)/eta_f);   % Fan's temperature ratio
    Tt1_3 = Tt2*tau_f;                                  % 
    Pt1_3 = Pt0*pi_d*pi_f;                              % 
    
    % 3.3. Secondary nozzle expansion
    Pt1_9 = delta_0*pi_d*pi_f*pi_ns*P0;
    M1_9 = sqrt((2/(gamma_c-1))*((Pt1_9/P0)^((gamma_c-1)/gamma_c)-1));
    Tt1_9 = Tt0*tau_f;
    
    if  M1_9 >= 1 % Chocked nozzle
        M1_9 = 1; 
        P1_9 = Pt1_9*(2/(gamma_c + 1))^(gamma_c/(gamma_c-1));
    else
        P1_9 = P0;
    end
    T1_9=Tt1_9 / (1 + ((gamma_c-1)/2)*M1_9^2);
    U1_9=sqrt(gamma_c*Rg*T1_9)*M1_9;
    
    % 3.3. Primary flow
    
    % 3.3.1. Low pressure compressor
    Pt2_5=Pt2 * pi_lpc;
    tau_lpc=1 + ((pi_lpc^((gamma_c - 1)/ gamma_c) - 1) / eta_lpc); %[Low pressure compressor's temperature ratio]
    Tt2_5=Tt2 * tau_lpc;
    
    % 3.3.2. High pressure compressor
    Pt3 = Pt2_5 * pi_hpc;
    tau_hpc = 1 + ((pi_hpc^((gamma_c - 1) / gamma_c) - 1) / eta_hpc); %[High pressure compressor's temperature ratio]
    Tt3 = Tt2_5 * tau_hpc;
    
    % 3.3.3. Burner
    f = (Cpt * Tt4 - Cpc * Tt3) / (eta_b * h - Cpt * Tt4); % Propellant fraction
    Pt4 = Pt3*pi_b;
    
    % 3.3.4. High pressure turbine
    tau_hpt = 1 - eta_mh^(-1) * (1+f)^(-1) * Cpc/Cpt * Tt2/Tt4 * tau_lpc * (tau_hpc - 1);%[High pressure turbine's temperature ratio]
    pi_hpt = (1+(tau_hpt-1)/eta_hpt)^(gamma_t/(gamma_t-1));%[High pressure turbine's pressure ratio]
    Pt4_5 = Pt4*pi_hpt;
    Tt4_5 = Tt4*tau_hpt;
    
    % 3.3.5. Low pressure turbine
    tau_lpt = 1-(1/(eta_ml*tau_hpt))*(Cpc/Cpt)*(Tt0/Tt4)*(1/(1+f))*((tau_lpc-1)+alpha*(tau_f-1));%[Low pressure turbine's temperature ratio]
    pi_lpt = (1+(tau_lpt-1)/eta_lpt)^(gamma_t/(gamma_t-1));%[High pressure turbine's pressure ratio]
    Pt5 = Pt4_5*pi_lpt;
    Tt5 = Tt4_5*tau_lpt;
    
    % 3.3.6. Primary nozzle expansion
    Tt9 = Tt5;          % Nozzle stagnation temperature [K]
    Pt9 = pi_np * Pt5;  % Nozzle stagnation pressure [Pa]
    M9 = sqrt(2/(gamma_t-1) * ((Pt9/P0)^((gamma_t-1)/gamma_t) - 1));    % Nozzle exhaust mach number [adim]
    P9 = 0;             % Nozzle exhaust pressure
    
    if  M9 >= 1 % Chocked nozzle
        M9 = 1; 
        P9 = Pt9*(2/(gamma_t + 1))^(gamma_t/(gamma_t-1));
    else
        P9 = P0;
    end
    T9 = Tt9/(1 + ((gamma_t-1)/2)*M9^2);    % Nozzle exhaust temperature [K]
    U9 = sqrt(gamma_t*Rg*T9)*M9;            % Nozzle exhaust velocity [m/s]


%% 4. Thrust for a TURBOFAN and final calcualtions

F_spec_p = (1+f)*U9 - U0 + ((1+f)*Rg*T9)/U9 * (1-P0/P9);
F_spec_s =  alpha*(U1_9 - U0) + alpha*(Rg*T1_9/U1_9)*(1- P0/P1_9);

F_spec_total = F_spec_p + F_spec_s;

%% 5. Results analysis

c_s = f/F_spec_total; % Specific fuel consumption
c_s_p = f/F_spec_p;
c_s_s = f/F_spec_s;

I_sp_p = (c_s_p * g)^(-1); % Spcific Impulse
I_sp_s = (c_s_s * g)^(-1); % Spcific Impulse
I_sp = (c_s * g)^(-1); % Spcific Impulse


% Outputs

% Temperature vector
T(1) = T0;
T(2) = Tt0;
T(3) = Tt1_3;
T(4) = Tt1_9;
T(5) = T1_9;
T(6) = Tt2;
T(7) = Tt2_5;
T(8) = Tt3;
T(9) = Tt4_5;
T(10) = Tt5; 
T(11) = Tt9;
T(12) = T9;

% Pressure vector
P(1) = Pt0;
P(2) = Pt1_3;
P(3) = Pt1_9;
P(4) = Pt2;
P(5) = Pt2_5;
P(6) = Pt3;
P(7) = Pt4;
P(8) = Pt4_5;
P(9) = Pt5;
P(10) = Pt9;

end


