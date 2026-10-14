clc; clear; close all;


%% Convergent-Divergent Nozzle - Setup
% --- Gas Properties (Air) ---
param.K = 1.2;                     % Specific heat ratio (gamma)
Mm = 0.024;                      % Molar mass [kg/mol]
R_univ = 8.3144;                   % Universal gas constant [J/(mol*K)]
param.R = R_univ / Mm;             % Specific gas constant for air [J/(kg*K)] -> ~287.05

% --- Reservoir / Stagnation Conditions ---
param.P0 = 40e+5;                % Total / stagnation pressure [Pa] 
param.T0 = 3000;                    % Total / stagnation temperature [K]

% --- Duct Geometry ---
param.L = 0.5;                     % Nozzle length [m]
param.xt = 0.15;                    % Throat position [m]
param.At = 0.0015;                  % Throat area [m^2]
param.A_1 = 0.02;                  % Inlet area [m^2]
param.A_2 = 0.025;                 % Exit area [m^2]

% Vectorized area definition as a Function Handle
param.A_func = @(x) (x <= param.xt) .* (param.At + (param.A_1 - param.At) * ((param.xt - x) / param.xt).^2) + ...
                    (x > param.xt) .* (param.At + (param.A_2 - param.At) * ((x - param.xt) / (param.L - param.xt)).^2);

% --- Nozzle operating limits (Calculated once to save computational time) ---
[M2_sub, M2_sup] = solve_Mach_from_Area(param.A_2/param.At, param.K);
param.M2_sub = M2_sub;
param.M2_sup = M2_sup;

% Division by the isentropic term (static pressure)
param.P2_sub = param.P0 / (1 + (param.K-1) / 2 * M2_sub^2)^(param.K / (param.K-1));
param.P2_sup = param.P0 / (1 + (param.K-1) / 2 * M2_sup^2)^(param.K / (param.K-1));

% Pressure for a normal shock located exactly at the exit section
param.P2_f = param.P2_sup * (1 + 2*param.K / (param.K+1) * (M2_sup^2 - 1));

% Factor for choked mass flow rate calculation
param.Gamma = (2 / (param.K + 1))^((param.K + 1) / (2 * (param.K - 1)));

%% ROCKET PARAMETERS
param.m_prop = 50;    % propellant mass [kg]
param.m_struct = 70;   % structure mass [kg]
param.m_pl = 15;       % payload mass [kg]
param.m0 = param.m_prop + param.m_pl + param.m_struct;
param.A = 0.05;          % section area [m^2]
param.gamma = 1.4;      % specific heat ratio of air

% Reference Mach vector and Cd
param.Mach_curve = [0, 0.8,  1.05, 1.2,  2.0,  3.0,  5.0,  10.0, 25.0];
param.Cd_curve   = [0.3, 0.3, 0.65, 0.55, 0.45, 0.35, 0.30, 0.30, 0.30];

% Calculate constant burning time based on initial mdot (since P0 is high, it's always choked in flight)
[initial_mdot, ~, ~] = evaluate_nozzle(101325, param);
param.tb = param.m_prop / initial_mdot;     % burning time [s]


%% SIMULATION
% initial conditions [altitude_0, velocity_0]
y0 = [0; 0];

% simulation interval time [t_start, t_end]
tspan = [0, 1000];

% Options to stop the simulation when the rocket touches the ground
options = odeset('Events', @ground_event);

% call to ode45 (function and parameters)
[t, y, te, ye, ~] = ode45(@(t, y) rocket_ode(t, y, param), tspan, y0, options);

% result extraction
h = y(:, 1); % first column: altitude
v = y(:, 2); % second column: velocity

% Display impact data on screen
fprintf('--- SIMULATION RESULTS ---\n');
if ~isempty(te)
    fprintf('Total flight time: %.2f s\n', te(1));
    fprintf('Maximum altitude: %.2f km\n', max(h)/1000);
    fprintf('Impact velocity: %.2f m/s (%.2f km/h)\n', ye(1,2), ye(1,2)*3.6);
else
    fprintf('No impact (Rocket at the ground).\n');
end


%% Post-Processing
P_history = zeros(length(t), 1);
T_history = zeros(length(t), 1);
regime_history = zeros(length(t), 1); 
x_shock_history = NaN(length(t), 1); 
q_history = zeros(length(t), 1);

for i = 1:length(t)
    [~, P_history(i), rho_i, ~] = standard_atmosphere(h(i));

    % Dynamic Pressure calculation: q = 0.5 * rho * v^2
    q_history(i) = 0.5 * rho_i * v(i)^2;

    if t(i) <= param.tb
        % Recalculate nozzle conditions for the specific historical altitude
        [mdot_i, V2_i, P2_i, regime_i, x_shock_i] = evaluate_nozzle(P_history(i), param);
        T_history(i) = mdot_i * V2_i + (P2_i - P_history(i)) * param.A_2;
        regime_history(i) = regime_i;
        x_shock_history(i) = x_shock_i;
    else
        T_history(i) = 0;
        regime_history(i) = 0;
    end
end

% Extraction of Max-Q point and corresponding conditions
[max_q, idx_max_q] = max(q_history);
t_max_q = t(idx_max_q);
h_max_q = h(idx_max_q);
v_max_q = v(idx_max_q);


%% Trajectory Graphs
figure;
subplot(2,1,1);
plot(t, h/1000, 'b', 'LineWidth', 1.5);
hold on;
if ~isempty(te)
    plot(te(1), ye(1,1)/1000, 'ro', 'MarkerFaceColor', 'r'); 
end
xlabel('Time [s]');
ylabel('Altitude [km]');
grid on;

subplot(2,1,2);
plot(t, v, 'r', 'LineWidth', 1.5);
hold on;
if ~isempty(te)
    plot(te(1), ye(1,2), 'ro', 'MarkerFaceColor', 'r'); 
end
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('Velocity Profile');
grid on;


%% Thrust and Pressure Graphs
figure;
subplot(2,1,1);
plot(t, T_history, 'b', 'LineWidth', 1.5);
title('Thrust Profile over Time');
xlabel('Time [s]'); 
ylabel('Thrust [N]');
grid on;

subplot(2,1,2);
plot(t, P_history, 'r', 'LineWidth', 1.5);
title('Ambient Pressure over Time');
xlabel('Time [s]'); 
ylabel('Pressure [Pa]');
grid on;


%% Nozzle Behaviour Graphs
% 1. Local pressure vs Physical threshold
figure;
subplot(3,1,1);
plot(t(t <= param.tb), P_history(t <= param.tb)/1e5, 'r', 'LineWidth', 2); hold on;
yline(param.P2_sub/1e5, '--w', 'P2_{sub} (Choking Limit)', 'LineWidth', 1.2);
yline(param.P2_f/1e5, '--b', 'P2_f (SW at the exit)', 'LineWidth', 1.2);
yline(param.P2_sup/1e5, '--g', 'P2_{sup} (Adapted Condition)', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Pressure [bar]');
title('Atmospheric Pressure Evolution and Nozzle Operative Threshold');
grid on;

% 2. Regime
subplot(3,1,2);
stairs(t(t <= param.tb), regime_history(t <= param.tb), 'b', 'LineWidth', 2);
yticks([1 2 3]);
yticklabels({'1: Subsonic', '2: Internal SW', '3: Supersonic'});
ylim([0.5 3.5]);
xlabel('Time [s]'); ylabel('Regime');
title('Nozzle State Transition during Ascending Phase');
grid on;

% 3. SW movement
subplot(3,1,3);
plot(t(t <= param.tb), x_shock_history(t <= param.tb), 'm', 'LineWidth', 2); hold on;
yline(param.xt, '--k', 'Throat (x_t)', 'LineWidth', 1.2);
yline(param.L, '--k', 'Exit (L)', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Shock Wave Position x [m]');
title('Shock Wave towards the exit');
ylim([param.xt*0.95, param.L*1.05]);
grid on;


%% Max-Q Graph
figure;
plot(t, q_history / 1000, 'k', 'LineWidth', 1.8); hold on;
plot(t_max_q, max_q / 1000, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

text(t_max_q + 1, max_q / 1000, sprintf('  Max-Q: %.2f kPa\n  t = %.1f s\n  h = %.1f km', ...
    max_q / 1000, t_max_q, h_max_q / 1000), 'FontSize', 10, 'FontWeight', 'bold');

xlabel('Time [s]');
ylabel('Dynamic Pressure q [kPa]');
title('Dynamic Pressure Profile and Max-Q Point');
grid on;


%% --- FUNCTIONS ---

function [mdot, V2, P2, regime, x_shock] = evaluate_nozzle(Ps, param)
    % Evaluates the nozzle regime based on external pressure (Ps)
    K = param.K; R = param.R; P0 = param.P0; T0 = param.T0; 
    At = param.At; A_2 = param.A_2; Gamma = param.Gamma;
    
    x_shock = NaN;

    if Ps >= param.P2_sub
        % 1. SUBSONIC REGIME
        regime = 1;
        M2 = sqrt(2 / (K-1) * ((P0 / Ps)^((K-1)/K) - 1));
        T2 = T0 / (1 + (K-1) / 2 * M2^2);
        a2 = sqrt(K * R * T2);
        A_cr = A_2 * M2 * ((2/(K+1)) * (1 + (K-1)/2 * (M2^2)))^(-(K+1)/(2*(K-1)));
        mdot = (A_cr * P0 / sqrt(R * T0)) * sqrt(K) * Gamma;
        P2 = Ps;

    elseif Ps <= param.P2_f
        % 2. SUPERSONIC REGIME (Choked with supersonic expansion)
        regime = 3;
        M2 = param.M2_sup;
        T2 = T0 / (1 + (K-1) / 2 * M2^2);
        a2 = sqrt(K * R * T2);
        mdot = (At * P0 / sqrt(R * T0)) * sqrt(K) * Gamma;
        P2 = param.P2_sup;
       
    else
        % 3. REGIME WITH NORMAL SHOCK WAVE IN THE DIVERGENT SECTION
        regime = 2;
        pos_min = param.xt;          
        pos_max = param.L;           
        tol = 10;              
        err = inf;             
        
        while err > tol
            pos = (pos_min + pos_max) / 2;
            A_local = param.A_func(pos);
            
            [~, Mm] = solve_Mach_from_Area(A_local/At, K);
            Mv = sqrt(((K-1)*Mm^2 + 2) / (2*K*Mm^2 - (K-1)));
            
            P0v = P0 * (((K+1)/(2*K*Mm^2 - (K-1)))^(1/(K-1))) * ((((K+1)*Mm^2)/((K-1)*Mm^2 + 2))^(K/(K-1)));
            
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

        x_shock = pos;
        T2 = T0 / (1 + (K-1) / 2 * M2^2);
        a2 = sqrt(K * R * T2);
        mdot = (At * P0 / sqrt(R * T0)) * sqrt(K) * Gamma;
    end
    
    V2 = M2 * a2; 
end

function dydt = rocket_ode(t, y, param)
    h = y(1); % y(1) = altitude h [m]
    v = y(2); % y(2) = velocity v [m/s]
    
    [~, P_amb, rho, a_local] = standard_atmosphere(h);

    R_Earth = 6371000; % radius of the Earth
    g0 = 9.80665; % gravity acceleration
    g = g0 * (R_Earth / (R_Earth + h))^2;

    if t <= param.tb
        % Evaluate dynamically the nozzle parameters at current external pressure
        [mdot, V2, P2] = evaluate_nozzle(P_amb, param);
        m = param.m0 - mdot * t;
        T = mdot * V2 + (P2 - P_amb) * param.A_2;
    else
        m = param.m_struct + param.m_pl; % Dry mass (constant)
        T = 0;
    end
    
    % instantaneous Mach number
    Mach = abs(v) / a_local;
    
    % linear interpolation of Cd based on Mach number
    Cd_actual = interp1(param.Mach_curve, param.Cd_curve, Mach, 'linear', 'extrap');
    
    % drag
    D = 0.5 * rho * v * abs(v) * param.A * Cd_actual;
    
    % instantaneous acceleration
    accel = (T - D) / m - g;       

    % derivatives of the system
    dydt = zeros(2,1);

    % ground constraint
    if h <= 0 && accel < 0
        dydt(1) = 0;
        dydt(2) = 0;
    else
        dydt(1) = v;                     % dh/dt = v
        dydt(2) = accel;                 % dv/dt = a
    end
end

function [value, isterminal, direction] = ground_event(~, y)
    value = y(1);      % check the altitude h = y(1)
    isterminal = 1;    % 1 = stop the simulation when the value is zero
    direction = -1;    % -1 = event only during the descent phase
end


function [Temp, P, rho, a_local] = standard_atmosphere(h)
    R = 287.05;
    gamma = 1.4;
    R_Earth = 6371000; % radius of the Earth
    g0 = 9.80665; % gravity acceleration
    g = g0 * (R_Earth / (R_Earth + h))^2;

    if h < 11000 % troposphere
        L = 0.0065;
        Temp = 288.15 - L * max(0, h);
        P = 101325 * (Temp / 288.15)^(g / (R * L));
        rho = P / (R * Temp);
    
    elseif h >= 11000 && h < 20000 % tropopause
        Temp = 216.65;
        P_11 = 22632.06;
        P = P_11 * exp(-(g / (R * Temp)) * (h - 11000));
        rho = P / (R * Temp);
    
    elseif h >= 20000 && h < 32000 % first layer stratosphere
        L = 0.001;
        Temp_20 = 216.65;
        P_20 = 5474.89;
        Temp = Temp_20 + L * (h - 20000);
        P = P_20 * (Temp / Temp_20)^(-(g / (R * L)));
        rho = P / (R * Temp);
    
    elseif h >= 32000 && h <= 47000 % second layer stratosphere
        L = 0.0028;
        Temp_32 = 228.65;
        P_32 = 868.02;
        Temp = Temp_32 + L * (h - 32000);
        P = P_32 * (Temp / Temp_32)^(-(g / (R * L)));
        rho = P / (R * Temp);
    
    else
        Temp = 273.15;
        rho = 0;
        P = 0;
    end

    a_local = sqrt(gamma * R * Temp);
end


function [M_sub, M_sup] = solve_Mach_from_Area(A_ratio, gamma)
    % A_ratio = A / A_star
    % Implicit equation definition g(M) = 0
    g = @(M) (1/M) * ((2/(gamma+1)) * (1 + (gamma-1)/2 * M^2))^((gamma+1)/(2*(gamma-1))) - A_ratio;
    
    % Subsonic root search
    M_sub = fzero(g, [1e-5, 1.0]);
    
    % Supersonic root search
    M_sup = fzero(g, [1.0, 10.0]);
end