clear; clc; close all
addpath(genpath(pwd));

params = config();

%% INITIAL STATE (same for both cases)
r0 = [params.R_e + 120e3; 0; 0];

gamma = deg2rad(-5);
v_mag = 7500;

r_hat = r0 / norm(r0);
ref = [0; 0; 1];
theta_hat = cross(ref, r_hat);
theta_hat = theta_hat / norm(theta_hat);

v0 = v_mag * (cos(gamma)*theta_hat + sin(gamma)*r_hat);

x0 = [r0; v0; 1; 0; zeros(3,1)];

nx = length(x0);
P0 = eye(nx);

%% =============================
% CASE 1: REACHABLE TARGET
% =============================
params.lat_target = deg2rad(0.3);
params.lon_target = deg2rad(9.2);

result_reach = simulate_mpc_ekf(x0, x0, P0, params);

r_reach = result_reach.x_true(:,1:3);
[lat_reach, lon_reach] = ecef_to_latlon(r_reach);

lat_f_r = lat_reach(end);
lon_f_r = lon_reach(end);

err_r = landing_error(lat_f_r, lon_f_r, ...
    params.lat_target, params.lon_target, params.R_e);

%% =============================
% CASE 2: UNREACHABLE TARGET
% =============================
params.lat_target = deg2rad(1.0);
params.lon_target = deg2rad(10.0);

result_unreach = simulate_mpc_ekf(x0, x0, P0, params);

r_unreach = result_unreach.x_true(:,1:3);
[lat_unreach, lon_unreach] = ecef_to_latlon(r_unreach);

lat_f_u = lat_unreach(end);
lon_f_u = lon_unreach(end);

err_u = landing_error(lat_f_u, lon_f_u, ...
    params.lat_target, params.lon_target, params.R_e);

%% =============================
% PRINT RESULTS
% =============================
fprintf('\n===== REACHABLE TARGET =====\n');
fprintf('Final Error: %.2f m\n', err_r);

fprintf('\n===== UNREACHABLE TARGET =====\n');
fprintf('Final Error: %.2f m\n', err_u);

%% =============================
% PLOT COMPARISON
% =============================
figure;

plot(rad2deg(lon_reach), rad2deg(lat_reach), 'b', 'LineWidth', 2); hold on;
plot(rad2deg(lon_unreach), rad2deg(lat_unreach), 'r--', 'LineWidth', 2);

% Plot targets
plot(9.2, 0.3, 'bo', 'MarkerFaceColor','b');
plot(10.0, 1.0, 'ro', 'MarkerFaceColor','r');

set(gca, 'FontSize', 24);

xlabel('Longitude (deg)', 'FontSize', 24);
ylabel('Latitude (deg)', 'FontSize', 24);
title('Reachable vs Unreachable Target Trajectories', 'FontSize', 24);

legend('Reachable Trajectory', ...
       'Unreachable Trajectory', ...
       'Reachable Target', ...
       'Unreachable Target','FontSize',24);

grid on;
axis equal;

exportgraphics(gcf, 'reachability.pdf', 'ContentType','vector');