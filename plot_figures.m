clear; clc; close all

% =============================
% GLOBAL FIGURE SETTINGS
% =============================
set(groot, 'defaultAxesFontSize', 24);
set(groot, 'defaultTextFontSize', 24);
set(groot, 'defaultLegendFontSize', 24);
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1);
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1);
set(groot, 'defaultLineLineWidth', 2);
% Create figures folder if it doesn't exist
fig_folder = fullfile(pwd, 'figures');
if ~exist(fig_folder, 'dir')
    mkdir(fig_folder);
end

%% LOAD DATA
addpath(genpath(pwd));
load('final_results_100MC.mat');   % loads mc, params, result_ekf

%% =============================
% EXTRACT NOMINAL DATA
% =============================
t = result_ekf.t;
x_true = result_ekf.x_true;
x_hat  = result_ekf.x_hat;
sigma  = result_ekf.sigma;

r_true = x_true(:,1:3);
v_true = x_true(:,4:6);

r_hat_est = x_hat(:,1:3);
v_hat_est = x_hat(:,4:6);

[lat_true, lon_true] = ecef_to_latlon(r_true);
[lat_hat,  lon_hat]  = ecef_to_latlon(r_hat_est);

alt   = vecnorm(r_true,2,2) - params.R_e;
speed = vecnorm(v_true,2,2);

%% =============================
% NOMINAL PLOTS
% =============================

% Altitude
fig1 = figure;
plot(t, alt/1000, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Altitude (km)');
title('Altitude vs Time');
grid on;
exportgraphics(fig1, fullfile(fig_folder, 'altitude.pdf'), 'ContentType','vector');

% Speed
fig2 = figure;
plot(t, speed, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Speed (m/s)');
title('Speed vs Time');
grid on;
exportgraphics(fig2, fullfile(fig_folder, 'speed.pdf'), 'ContentType','vector');
% Ground Track
fig3 = figure;
plot(rad2deg(lon_true), rad2deg(lat_true), 'b', 'LineWidth', 2); hold on;
plot(rad2deg(lon_hat),  rad2deg(lat_hat),  'r--', 'LineWidth', 2);
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('Ground Track (True vs Estimated)');
legend('True','Estimated');
grid on;
exportgraphics(fig3, fullfile(fig_folder, 'ground_track.pdf'), 'ContentType','vector');

% Bank Angle
t_sigma = t(1:length(sigma));
fig4 = figure;
plot(t_sigma, rad2deg(sigma), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Bank Angle (deg)');
title('Bank Angle History');
grid on;
exportgraphics(fig4, fullfile(fig_folder, 'bank_angle.pdf'), 'ContentType','vector');

%% =============================
% NOMINAL FINAL ERROR
% =============================
lat_f = lat_true(end);
lon_f = lon_true(end);

[downrange, crossrange] = compute_range_errors( ...
    r_true(end,:)', params.lat_target, params.lon_target, params.R_e);

total_error = landing_error(lat_f, lon_f, ...
    params.lat_target, params.lon_target, params.R_e);

fprintf('\n===== NOMINAL RESULT =====\n');
fprintf('Total Error: %.2f m\n', total_error);

%% =============================
% MONTE CARLO STATS
% =============================
fprintf('\n===== MONTE CARLO STATS =====\n');

fprintf('Mean Error      : %.2f m\n', mean(mc.total_error));
fprintf('Std Dev Error   : %.2f m\n', std(mc.total_error));
fprintf('Median Error    : %.2f m\n', median(mc.total_error));
fprintf('Max Error       : %.2f m\n', max(mc.total_error));
fprintf('Min Error       : %.2f m\n', min(mc.total_error));

fprintf('Mean Downrange  : %.2f m\n', mean(mc.downrange));
fprintf('Mean Crossrange : %.2f m\n', mean(mc.crossrange));

fprintf('\n===== 3-SIGMA =====\n');
fprintf('3σ Downrange  : %.2f m\n', 3*std(mc.downrange));
fprintf('3σ Crossrange : %.2f m\n', 3*std(mc.crossrange));
fprintf('3σ Total Err  : %.2f m\n', 3*std(mc.total_error));

%% =============================
% CONSTRAINT VIOLATIONS
% =============================
p_q    = mean(mc.q_peak > params.q_max);
p_g    = mean(mc.g_peak > params.g_max);
p_heat = mean(mc.heat_peak > params.heat_max);

fprintf('\n===== CONSTRAINT VIOLATIONS =====\n');
fprintf('P(q > q_max)       : %.3f\n', p_q);
fprintf('P(g > g_max)       : %.3f\n', p_g);
fprintf('P(heat > heat_max) : %.3f\n', p_heat);

%% =============================
% MONTE CARLO PLOTS
% =============================

% Histogram
fig5 = figure;
histogram(mc.total_error/1000, 15);
xlabel('Landing Error (km)');
ylabel('Frequency');
title('Landing Error Distribution');
grid on;
exportgraphics(fig5, fullfile(fig_folder, 'mc_histogram.pdf'), 'ContentType','vector');

% Dispersion
fig6 = figure;
scatter(mc.downrange/1000, mc.crossrange/1000, 30, 'filled');
xlabel('Downrange Error (km)');
ylabel('Crossrange Error (km)');
title('Landing Dispersion');
grid on;
axis equal;
exportgraphics(fig6, fullfile(fig_folder, 'mc_dispersion.pdf'), 'ContentType','vector');

% 3σ Ellipse
fig7 = figure;
scatter(mc.downrange/1000, mc.crossrange/1000, 30, 'filled'); hold on;

mu = [mean(mc.downrange); mean(mc.crossrange)];
Sigma = cov(mc.downrange, mc.crossrange);

theta = linspace(0, 2*pi, 200);
circle = [cos(theta); sin(theta)];

[U,S,~] = svd(Sigma);
ellipse = U * sqrt(S) * circle * 3;

plot((ellipse(1,:) + mu(1))/1000, (ellipse(2,:) + mu(2))/1000, ...
     'r', 'LineWidth', 2);

xlabel('Downrange Error (km)');
ylabel('Crossrange Error (km)');
title('Landing Dispersion with 3σ Ellipse');
legend('Samples','3σ ellipse');
grid on;
axis equal;

exportgraphics(fig7, fullfile(fig_folder, 'mc_ellipse.pdf'), 'ContentType','vector');

% CDF
fig8 = figure;
[f, x] = ecdf(mc.total_error);
plot(x/1000, f, 'LineWidth', 2);
xlabel('Landing Error (km)');
ylabel('Cumulative Probability');
title('CDF of Landing Error');
grid on;
exportgraphics(fig8, fullfile(fig_folder, 'mc_cdf.pdf'), 'ContentType','vector');

set(gcf, 'Renderer', 'painters');