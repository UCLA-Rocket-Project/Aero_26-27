clear; clc;

%INPUTS
%--------------------------------------------------------------------------------------
%NC
L_nc_in = 30;                      % Length of the nose cone [in]
L_shldr_in = 3;                    % [in]
layers = 5;                        % #fiber glass layers
D_nc = 8;                          % Nose cone base diameter [in] 
wall_thickness = layers*0.0004;    % Wall thickness of the nose cone [m]
rho_fiberglass = 1900;             % Density of fiberglass [kg/m³] applys to 9oz s glass, confirm if this is backed up with real world or just online numbers
rho_alum = 2810;                   % 7075 aluminum density [kg/m^3]
%SPIKE
L_spk_in = 8;                      % [in]
D_spk = 0.75;                      % [in]
%TOTAL
a = 10;                            % rocket acceleration [m/s^2]
g = 9.81;                          % graviational acceleration[m/s^2]
%--------------------------------------------------------------------------------------

%convert
R_nc = (D_nc/2)*0.0254;        %[m]
L_nc = L_nc_in*0.0254;         %[m]
L_shldr = L_shldr_in*0.0254;   %[m]
L_spk = L_spk_in*0.0254;       %[m]
R_spk = D_spk*0.0254/2;        %[m]

% CALCULATE NOSE CONE MASS
% Von Karman profile: y(x) = (R/sqrt(pi)) * sqrt(acos(1-2*x/L) - sin(2*acos(1-2*x/L))/2
% We need to integrate to find volume

% Discretize nose cone into slices
n_slices = 10000;
x = linspace(0, L_nc, n_slices);

% Calculate radius at each position using Von Karman equation
theta = acos(1 - (2*x)/L_nc);
y = (R_nc/sqrt(pi)) * sqrt(theta - sin(2*theta)/2);

% Calculate radi
% Shell volume = Outer - Inner
r_outer = y;
r_inner = max(0, y - wall_thickness); % Don't go negative near tip

% Total nc volume
V_ncone = trapz(x, pi * (r_outer.^2 - r_inner.^2));
V_shldr = pi*((R_nc^2)-((R_nc-wall_thickness)^2))*L_shldr;
V_nc = V_shldr +V_ncone;

% Total nc mass
mass_nc = V_nc * rho_fiberglass; % [kg]

% CALCULATE SPIKE MASS
% volume spike
V_shaft = pi*(R_spk^2)*L_spk;
V_tip = (pi*2/3)*(R_spk^3);
%V_connection =
%V_anchor =
V_spk = V_shaft+V_tip;           %[m^3]

% mass spike
mass_spk = V_spk * rho_alum;     % [kg]

%total upward applied force
mass_total = mass_spk+mass_nc;   % [kg]
Fa_spk = mass_spk*(a+g);         % [N]
Fa_nc = mass_nc*(a+g);           % [N]
Fa_total = mass_total*(a+g);     % [N]

% print values
fprintf("Spike mass: %.3f kg(%.3f lbm) *not including anchor/connector\n", mass_spk, mass_spk*2.20462);
fprintf("Nose cone mass: %.3f kg (%.3f lbm) *not including thickened\n", mass_nc, mass_nc*2.20462);
fprintf("Total mass: %.3f kg (%.3f lbm)\n", mass_total, mass_total*2.20462);
fprintf("Aplied upward force:%.3f N (%.3f pound-force)\n", Fa_total, Fa_total*0.224809);
