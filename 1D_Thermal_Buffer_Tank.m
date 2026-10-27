% 1D Thermal Buffer Tank Simulation with Baffle Effect
% Compares 25s transient numerical results against STAR-CCM+ export
% Author Ameya Desai

clear; clc; close all;

%% 1. Physical and Domain Parameters
H = 1.5;                  % Tank height (m)
N = 200;                  % Number of spatial nodes
dz = H / (N - 1);         % Spatial step size (m)
z = linspace(0, H, N)';   % Spatial grid (m)

t_span = [0, 25.0];       % Simulation time (s)

% Fluid Properties (Water approximation)
alpha_base = 1.4e-7;      % Base thermal diffusivity (m^2/s)
v_inlet = 0.05;           % Advection velocity matched to CFD front

% Temperatures (K) 
T_hot = 333.0;            % Inlet fluid temperature (Top)
T_cold = 293.0;           % Initial tank temperature

%% 2. Baffle Modeling (Effective Diffusivity Modification)
baffle_height = 0.5;      % Baffle location (m)
alpha_eff = alpha_base * ones(N, 1);

for i = 1:N
    if abs(z(i) - baffle_height) < 0.05
        alpha_eff(i) = alpha_base * 0.1; % Suppressed mixing zone at baffle
    end
end

%% 3. Initial Conditions
T_initial = T_cold * ones(N, 1);
T_initial(N) = T_hot;     % Enforce Dirichlet BC in the initial state vector

%% 4. Solve 1D Advection-Diffusion using ode45
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
[t, T_out] = ode45(@(t, T) tank_pde(t, T, dz, v_inlet, alpha_eff), t_span, T_initial, options);

T_matlab_final = T_out(end, :)';

%% 5. Plot Comparison (Report-Formatted)
figure('Name', 'Thermal Stratification Comparison', 'Color', 'w');
hold on; grid on;

% Apply white background and black axes for report compatibility
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 11);

% Plot MATLAB Data
plot(T_matlab_final, z, 'b-', 'LineWidth', 2, 'DisplayName', 'MATLAB (1D PDE)');

% Attempt to load and plot STAR-CCM+ Data
try
    star_data = readmatrix('temp_starccm.csv'); 
    z_star = star_data(:, 1);
    T_star = star_data(:, 2);
    
    plot(T_star, z_star, 'r-', 'LineWidth', 2, 'DisplayName', 'STAR-CCM+');
catch
    warning('temperature_profile.csv not found. Plotting only MATLAB results.');
end

% Formatting Labels and Titles
xlabel('Temperature (K)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
ylabel('Tank Height, Z (m)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
title('Thermal Buffer Tank Profile at t = 25.0 s', 'FontSize', 14, 'Color', 'k');

% Formatting Legend
lgd = legend('Location', 'best', 'FontSize', 11);
set(lgd, 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');

%% Local Function for PDE (Method of Lines)
function dTdt = tank_pde(~, T, dz, v, alpha_eff)
    N = length(T);
    dTdt = zeros(N, 1);
    
    for i = 2:N-1
        % Top-down flow upwind scheme
        advection = v * (T(i) - T(i+1)) / dz; 
        diffusion = alpha_eff(i) * (T(i+1) - 2*T(i) + T(i-1)) / (dz^2);
        dTdt(i) = -advection + diffusion;
    end
    
    % Boundary Conditions
    dTdt(N) = 0;                       % Top boundary remains at T_hot
    dTdt(1) = -v * (T(1) - T(2)) / dz; % Bottom boundary outflow
end