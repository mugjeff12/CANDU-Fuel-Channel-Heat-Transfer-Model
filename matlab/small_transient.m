%% 02_small_transient.m  —  CORRECTED
% Small transient via Method of Lines for a CANDU-like fuel channel.
% Adds the missing thermal capacitance term (rho*A*dx*cp) so temperatures
% evolve physically (no blow-up). Uses upwind advection and a step in q'.
%
% Save as: matlab/02_small_transient.m

clear; clc;

%% ---------------- Geometry & constants ----------------
L   = 6.0;          % m   channel length
Dh  = 0.012;        % m   hydraulic diameter
Pch = 10e6;         % Pa  pressure
Tin = 560.0;        % K   inlet coolant temperature
eps = 0.0;          % m   roughness (unused here; transient is thermal only)

N  = 60;                               % axial cells
x  = linspace(0, L, N).';              % cell centers
dx = x(2) - x(1);

Pwet  = pi*Dh;                          % m  wetted perimeter
Aflow = pi*(Dh/2)^2;                    % m^2 flow area

%% ---------------- Operating point & transient ----------------
mdot_base = 2.0;         % kg/s  base flow
q_base    = 15.5e3;      % W/m   base linear heat rate (≈ +10 K margin case)

% Step in power at t = t_step  (set q_step = 1.0 for no step)
q_step = 1.6;            % multiplier after the step (e.g., 1.6 => +60% power)
t_step = 10.0;           % s

% If you prefer a flow step instead, set:
use_flow_step = false;   % true => step mdot instead of q'
mdot_step = 1.4;         % multiplier when use_flow_step = true

%% ---------------- Initial condition ----------------
% Start close to steady by using a mild axial gradient (optional).
Tc0 = Tin * ones(N,1);   % K (simple and robust)

%% ---------------- Integrate in time ----------------
tspan = [0 40];          % seconds
opts  = odeset('RelTol',1e-6,'AbsTol',1e-7);
[t_sol, Tc_sol] = ode15s(@rhs, tspan, Tc0, opts);

%% ---------------- Post-process & plots ----------------
Tout = Tc_sol(:, end);

% 1) Outlet response
figure('Color','w');
plot(t_sol, Tout, 'LineWidth', 2);
grid on; xlabel('Time (s)'); ylabel('Outlet coolant T_{out} (K)');
title('Transient outlet response to step');

% 2) Space–time map
figure('Color','w');
imagesc(x, t_sol, Tc_sol); set(gca,'YDir','normal');
colorbar; xlabel('x (m)'); ylabel('time (s)');
title('Coolant temperature T_c(x,t)');

%% ---------------- RHS: dT/dt = (conv + source)/cap ----------------
function dTdt = rhs(t, Tc)
    % pull parameters from base workspace (simplifies a single-file script)
    L   = evalin('base','L'); Dh  = evalin('base','Dh');  Pch = evalin('base','Pch');
    Tin = evalin('base','Tin'); N   = evalin('base','N');  dx  = evalin('base','dx');
    Pwet= evalin('base','Pwet'); Aflow = evalin('base','Aflow');
    mdot_base = evalin('base','mdot_base'); q_base = evalin('base','q_base');
    q_step = evalin('base','q_step'); t_step = evalin('base','t_step');
    use_flow_step = evalin('base','use_flow_step'); mdot_step = evalin('base','mdot_step');

    % step logic
    if t < t_step
        q_lin = q_base;
        mdot  = mdot_base;
    else
        q_lin = q_base * (use_flow_step*1 + (~use_flow_step)*q_step) + ...
                q_base * (use_flow_step*0);
        mdot  = mdot_base * (use_flow_step*mdot_step + (~use_flow_step)*1);
    end

    % Heat input per unit area
    qpp = q_lin / Pwet;   % W/m^2

    % Film-T property update and local h (carry downstream)
    rho = zeros(N,1); cp = zeros(N,1);
    h_prev = [];
    for i = 1:N
        if isempty(h_prev)
            pr0 = water_props_T(Tc(1), Pch);
            [h_prev, ~, ~, ~] = h_from_props(pr0.rho, pr0.mu, pr0.k, pr0.cp, Dh, mdot, Aflow);
        end
        Tw_est = Tc(i) + qpp/max(h_prev, 1e-3);
        Tf  = 0.5*(Tc(i) + Tw_est);
        pr  = water_props_T(Tf, Pch);
        [h_i, ~, ~, ~] = h_from_props(pr.rho, pr.mu, pr.k, pr.cp, Dh, mdot, Aflow);
        rho(i) = pr.rho; cp(i) = pr.cp;
        h_prev = h_i; % march downstream
    end

    % Upwind first-order convection (Dirichlet inlet Tin)
    dTdx = zeros(N,1);
    for i = 1:N
        if i == 1
            dTdx(i) = (Tc(i) - Tin)/dx;
        else
            dTdx(i) = (Tc(i) - Tc(i-1))/dx;
        end
    end

    % --- Correct thermal capacity form ---
    % cap [J/K per cell] = rho * A * dx * cp
    cap = rho .* Aflow .* dx .* cp;

    % convective energy flux gradient [W] = mdot*cp * dT/dx
    conv_W = (mdot .* cp) .* dTdx;

    % linear heat input into each cell [W] = q'' * P * dx
    source_W = (qpp .* Pwet .* dx) .* ones(N,1);

    % time derivative [K/s]
    dTdt = -(conv_W - source_W) ./ cap;  % minus sign because dTdx uses upwind
end

%% ---------------- Helpers (same trends as Python props) ----------------
function [h, Re, Nu, v] = h_from_props(rho, mu, k, cp, Dh, mdot, Aflow)
    v  = mdot/(rho*Aflow);
    Re = rho*v*Dh/mu;
    Pr = cp*mu/k;
    Nu = 0.023*(Re^0.8)*(Pr^0.4);
    h  = Nu*k/Dh;
end

function props = water_props_T(T, P)
% Educational water-like property trends near 10 MPa (illustrative).
    rho = 700 - 0.2*(T-560); rho = max(rho,200);
    mu  = 1.2e-4 * exp(-0.012*(T-560));    % Pa·s
    k   = 0.65 - 1.0e-4*(T-560);           % W/m-K
    cp  = 5000 + 0.5*(T-560);              % J/kg-K
    props = struct('rho',rho,'mu',mu,'k',k,'cp',cp);
end
