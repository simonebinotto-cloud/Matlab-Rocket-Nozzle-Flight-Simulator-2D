clc; clear; close all;

% --- Gas Properties (Air) ---
K = 1.4;                     % Specific heat ratio (gamma)
Mm = 0.02897;                % Molar mass [kg/mol]
R_univ = 8.3144;             % Universal gas constant [J/(mol*K)]
R = R_univ / Mm;             % Specific gas constant for air [J/(kg*K)] -> ~287.05

% --- Reservoir / Stagnation Conditions ---
P0 = 1013250;                % Total / stagnation pressure [Pa] (10 bar)
T0 = 300;                    % Total / stagnation temperature [K]
Ps = 0.60 * P0;              % External Pressure

% --- Duct Geometry ---
L = 1.0;                     % Nozzle length [m]
xt = 0.4;                    % Throat position [m]
At = 0.01;                   % Throat area [m^2]
A_1 = 0.02;                  % Inlet area [m^2]
A_2 = 0.015;                 % Exit area [m^2]

% Vectorized area definition as a Function Handle
A_func = @(x) (x <= xt) .* (At + (A_1 - At) * ((xt - x) / xt).^2) + ...
              (x > xt) .* (At + (A_2 - At) * ((x - xt) / (L - xt)).^2);

% --- Nozzle operating limits ---
[M2_sub, M2_sup] = solve_Mach_from_Area(A_2/At, K);

% Division by the isentropic term (static pressure)
P2_sub = P0 / (1 + (K-1) / 2 * M2_sub^2)^(K / (K-1));
P2_sup = P0 / (1 + (K-1) / 2 * M2_sup^2)^(K / (K-1));

% Pressure for a normal shock located exactly at the exit section
P2_f = P2_sup * (1 + 2*K / (K+1) * (M2_sup^2 - 1));

% Vandenkerckhove factor for choked mass flow rate calculation
Gamma = (2 / (K + 1))^((K + 1) / (2 * (K - 1)));

% --- Flow Regimes Analysis ---
if Ps >= P2_sub
    % 1. SUBSONIC REGIME
    disp('Flow Regime: Subsonic (Unchoked nozzle)');
    
    M2 = sqrt(2 / (K-1) * ((P0 / Ps)^((K-1)/K) - 1));
    T2 = T0 / (1 + (K-1) / 2 * M2^2);
    a2 = sqrt(K * R * T2);
    A_cr = A_2 * M2 * ((2/(K+1)) * (1 + (K-1)/2 * (M2^2)))^(-(K+1)/(2*(K-1)));
    mdot = (A_cr * P0 / sqrt(R * T0)) * sqrt(K) * Gamma;
    P2 = Ps;

elseif Ps <= P2_f
    % 2. SUPERSONIC REGIME (Choked with supersonic expansion)
    disp('Flow Regime: Supersonic / Under-expanded');
    M2 = M2_sup;
    T2 = T0 / (1 + (K-1) / 2 * M2^2);
    a2 = sqrt(K * R * T2);
    
    % Mass flow rate is maximized and determined by the throat
    mdot = (At * P0 / sqrt(R * T0)) * sqrt(K) * Gamma;
    P2 = P2_sup;
   
elseif Ps > P2_f && Ps < P2_sub
    % 3. REGIME WITH NORMAL SHOCK WAVE IN THE DIVERGENT SECTION
    disp('Flow Regime: Normal shock wave in the divergent section');
    
    % Bisection method settings
    pos_min = xt;          
    pos_max = L;           
    tol = 10;              
    err = inf;             
    
    while err > tol
        pos = (pos_min + pos_max) / 2;
        
        % Retrieve local area using the function handle
        A_local = A_func(pos);
        
        [~, Mm] = solve_Mach_from_Area(A_local/At, K);
        Mv = sqrt(((K-1)*Mm^2 + 2) / (2*K*Mm^2 - (K-1)));
        
        P0v = P0 * (((K+1)/(2*K*Mm^2 - (K-1)))^(1/(K-1))) * ((((K+1)*Mm^2)/((K-1)*Mm^2 + 2))^(K/(K-1)));
        Pv = P0v / ((1 + (K-1)/2 * Mv^2)^(K/(K-1)));
        
        A_cr = A_local * Mv * ((2/(K+1)) * (1 + (K-1)/2 * Mv^2))^(-(K+1)/(2*(K-1)));
        [M2, ~] = solve_Mach_from_Area(A_2/A_cr, K);
        
        P2 = P0v / ((1 + (K-1)/2 * M2^2)^(K/(K-1)));
        err = abs(P2 - Ps);
        
        if P2 < Ps
            pos_max = pos;
        else
            pos_min = pos;
        end
    end
    
    T2 = T0 / (1 + (K-1) / 2 * M2^2);
    a2 = sqrt(K * R * T2);
    
    % Since the throat is choked, mass flow rate is at its maximum
    mdot = (At * P0 / sqrt(R * T0)) * sqrt(K) * Gamma;
    fprintf(' -> Shock wave position found at x = %.4f m\n', pos);
end

% --- THRUST CALCULATION ---
V2 = M2 * a2;                                % Exit velocity [m/s]
Thrust = mdot * V2 + (P2 - Ps) * A_2;        % Total thrust [N]

fprintf(' -> Mass Flow Rate = %.3f kg/s\n', mdot);
fprintf(' -> Exit Velocity = %.2f m/s\n', V2);
fprintf(' -> Total Thrust = %.2f N\n', Thrust);


function [M_sub, M_sup] = solve_Mach_from_Area(A_ratio, gamma)
    % A_ratio = A / A_star
    
    % Implicit equation definition g(M) = 0
    g = @(M) (1/M) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * M^2))^((gamma+1)/(2*(gamma-1))) - A_ratio;
    
    % Subsonic root search
    M_sub = fzero(g, [1e-5, 1.0]);
    
    % Supersonic root search
    M_sup = fzero(g, [1.0, 10.0]);
end