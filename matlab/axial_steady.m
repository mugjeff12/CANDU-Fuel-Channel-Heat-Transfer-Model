%% 01_axial_steady.m
% CANDU Fuel Channel – steady 1-D model with film-T properties,
% Dittus–Boelter h, Churchill friction, and margin to boiling.
% Save this file as: matlab/01_axial_steady.m  (R2016b+ for local functions)

clear; clc;

%% ---------------- Base inputs ----------------
L     = 6.0;         % m  channel length
Dh    = 0.012;       % m  hydraulic diameter
Tin   = 560.0;       % K  inlet coolant temperature
P     = 10e6;        % Pa pressure
mdot  = 2.0;         % kg/s mass flow rate
q_lin = 15.5e3;      % W/m linear heat rate (e.g., tuned for +10 K margin)
eps   = 0.0;         % m absolute roughness (0 ~ smooth metal)

N  = 200;                        % axial nodes
x  = linspace(0, L, N+1).';      % column vector
dx = x(2)-x(1);

Pwet  = pi*Dh;
Aflow = pi*(Dh/2)^2;
qpp   = q_lin / Pwet;            % W/m^2

%% ---------------- Arrays ----------------
Tc = zeros(N+1,1); Tc(1) = Tin;
Tw = zeros(N+1,1);
h  = zeros(N+1,1);
Re = zeros(N+1,1);
Nu = zeros(N+1,1);
ff = zeros(N+1,1);
dP = zeros(N+1,1);               % cumulative Pa

%% ---------------- Inlet film estimate ----------------
props0 = water_props_T(Tin, P);
[h0, Re0, Nu0, v0] = h_from_props(props0.rho, props0.mu, props0.k, props0.cp, Dh, mdot, Aflow);
Tw(1) = Tc(1) + qpp/h0;
h(1)  = h0;  Re(1)=Re0;  Nu(1)=Nu0;
ff(1) = f_churchill(Re(1), Dh, eps);

%% ---------------- Axial march with film-T iteration ----------------
for i = 1:N
    % 1) Advance coolant using cp at current film estimate
    props_i = water_props_T(0.5*(Tc(i)+Tw(i)), P);
    Tc(i+1) = Tc(i) + (q_lin/(mdot*props_i.cp))*dx;

    % 2) Film-T iteration (few fixed-point passes)
    Tc_loc = Tc(i+1);
    Tw_loc = Tc_loc + qpp/max(h(i),1e-3); % start from previous h
    pr = props_i;
    for k = 1:3
        Tf = 0.5*(Tc_loc + Tw_loc);
        pr = water_props_T(Tf, P);
        [h_loc, Re_loc, Nu_loc, v_loc] = h_from_props(pr.rho, pr.mu, pr.k, pr.cp, Dh, mdot, Aflow);
        Tw_loc = Tc_loc + qpp/h_loc;
    end

    Tw(i+1) = Tw_loc;  h(i+1) = h_loc;  Re(i+1) = Re_loc;  Nu(i+1) = Nu_loc;
    ff(i+1) = f_churchill(Re_loc, Dh, eps);

    % 3) Hydraulics
    dPdx = ff(i+1)*(pr.rho*v_loc^2/2)/Dh;
    dP(i+1) = dP(i) + dPdx*dx;
end

%% ---------------- Margin to boiling ----------------
Tsat = tsat_water(P);                 % K (CoolProp if present; else ~584 K)
DeltaT_sat = Tsat - Tw;
[min_margin, ix_min] = min(DeltaT_sat);
x_min = x(ix_min);

%% ---------------- Save CSV (data/axial_steady_base.csv) ----------------
base    = fileparts(fileparts(pwd));  % repo root if you're in matlab/
dataDir = fullfile(base,'data'); if ~exist(dataDir,'dir'), mkdir(dataDir); end
T = table(x, Tc, Tw, h, Re, Nu, ff, dP, repmat(Tsat,size(x)), DeltaT_sat, ...
    'VariableNames', {'x_m','T_c_K','T_w_K','h_W_m2K','Re','Nu','f','dP_cum_Pa','T_sat_K','DeltaT_sat_K'});
writetable(T, fullfile(dataDir,'axial_steady_base.csv'));

%% ---------------- Plots (figures/*.png) ----------------
figDir = fullfile(base,'figures'); if ~exist(figDir,'dir'), mkdir(figDir); end

% 1) Coolant, wall, Tsat
f1 = figure('Color','w');
plot(x,Tc,'LineWidth',2); hold on;
plot(x,Tw,'--','LineWidth',2);
plot(x, Tsat*ones(size(x)),'-.','LineWidth',1.5);
grid on; xlabel('Axial position x (m)'); ylabel('Temperature (K)');
title('Coolant, Wall, and Saturation Temperature vs x');
legend('Coolant T_c','Wall T_w','T_{sat} (≈)','Location','best');
exportgraphics(f1, fullfile(figDir,'Tx_coolant_wall_Tsat.png'), 'Resolution',200);

% 2) ΔP
f2 = figure('Color','w');
plot(x, dP/1e6, 'LineWidth',2);
grid on; xlabel('Axial position x (m)'); ylabel('\DeltaP (MPa)');
title('Cumulative Pressure Drop vs x (Churchill friction)');
exportgraphics(f2, fullfile(figDir,'dP_x_filmT.png'), 'Resolution',200);

% 3) Margin + annotation
f3 = figure('Color','w');
plot(x, DeltaT_sat, 'LineWidth',2); hold on;
yline(0,'--'); scatter(x_min, min_margin, 40, 'filled');
text(x_min, min_margin, sprintf('  min %.1f K @ x=%.2f m',min_margin,x_min), 'VerticalAlignment','bottom');
grid on; xlabel('Axial position x (m)'); ylabel('\DeltaT_{sat} (K)');
title('Margin to Boiling vs x');
exportgraphics(f3, fullfile(figDir,'margin_to_boiling.png'), 'Resolution',200);

%% ---------------- Console summary ----------------
fprintf('=== Axial Steady (MATLAB) – Film-T + Churchill + Margin ===\n');
fprintf('T_out = %.2f K, T_w,out = %.2f K, T_sat ≈ %.2f K\n', Tc(end), Tw(end), Tsat);
fprintf('h_in = %.0f W/m^2-K, h_out = %.0f W/m^2-K\n', h(1), h(end));
fprintf('Re range ~ %.2e → %.2e\n', Re(1), Re(end));
fprintf('Friction: Churchill, eps = %g m\n', eps);
fprintf('Total ΔP = %.3f MPa over L = %.1f m\n', dP(end)/1e6, L);
fprintf('Min ΔT_sat = %.1f K at x = %.2f m\n', min_margin, x_min);
fprintf('CSV:    %s\n', fullfile(dataDir,'axial_steady_base.csv'));
fprintf('Figures:%s\n        %s\n        %s\n', ...
    fullfile(figDir,'Tx_coolant_wall_Tsat.png'), ...
    fullfile(figDir,'dP_x_filmT.png'), ...
    fullfile(figDir,'margin_to_boiling.png'));

%% ---------------- Local functions ----------------
function [h, Re, Nu, v] = h_from_props(rho, mu, k, cp, Dh, mdot, Aflow)
    v  = mdot/(rho*Aflow);
    Re = rho*v*Dh/mu;
    Pr = cp*mu/k;
    Nu = 0.023*(Re^0.8)*(Pr^0.4);
    h  = Nu*k/Dh;
end

function f = f_churchill(Re, Dh, eps)
    if Re < 1e-12, f = 0.0; return; end
    A = (2.457*log((7.0/Re)^0.9 + 0.27*(eps/Dh)))^16;
    B = (37530.0/Re)^16;
    f = 8.0 * ((8.0/Re)^12 + 1.0/((A+B)^1.5))^(1/12);
end

function Tsat = tsat_water(P)
% Tries CoolProp if available; else returns ~584 K at 10 MPa (illustrative).
    Tsat = 584.0;
    try
        if exist('CoolProp.PropsSI','file')
            Tsat = CoolProp.PropsSI('T','P',P,'Q',0,'Water');
        end
    catch
        % keep placeholder
    end
end

function props = water_props_T(T, P)
% Educational property fit: water-like at high pressure.
% Returns struct with rho [kg/m3], mu [Pa.s], k [W/m-K], cp [J/kg-K]
% Simple temperature trends around 550–650 K range (not D2O, illustrative).
    rho = 700 - 0.2*(T-560);         % very rough trend
    rho = max(rho, 200);
    mu  = 1.2e-4 * exp(-0.012*(T-560));  % Pa.s, decreasing with T
    k   = 0.65 - 1.0e-4*(T-560);     % W/m-K
    cp  = 5000 + 0.5*(T-560);        % J/kg-K
    props = struct('rho',rho,'mu',mu,'k',k,'cp',cp);
end
