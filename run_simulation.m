clear; clc; close all
addpath(genpath(pwd));

params = config();

% Initial true state
r0 = [params.R_e + 120e3; 0; 0];

gamma = deg2rad(-5);
v_mag = 7500;

r_hat = r0 / norm(r0);
ref = [0; 0; 1];
theta_hat = cross(ref, r_hat);
theta_hat = theta_hat / norm(theta_hat);

v0 = v_mag * (cos(gamma)*theta_hat + sin(gamma)*r_hat);

beta_scale = 1;
rho_bias   = 0;
wind       = [0; 0; 0];

x0_true = [r0; v0; beta_scale; rho_bias; wind];

% Initial estimate
x0_hat = x0_true;
x0_hat(1:3) = x0_hat(1:3) + 1000*randn(3,1);
x0_hat(4:6) = x0_hat(4:6) + 10*randn(3,1);

nx = length(x0_true);
P0 = diag([ ...
    1e6 * ones(3,1); ...
    1e2 * ones(3,1); ...
    0.1; ...
    0.1; ...
    100 * ones(3,1) ...
]);

% Run one EKF-MPC simulation
result_ekf = simulate_mpc_ekf(x0_true, x0_hat, P0, params);

% Monte Carlo
Nsim =100;
mc = monte_carlo_ekf(x0_true, params, Nsim);



%% =============================
% EXTRACT STATES
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

figure;
plot(t, alt/1000, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Altitude (km)');
title('Altitude vs Time');
grid on;

figure;
plot(t, speed, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Speed (m/s)');
title('Speed vs Time');
grid on;


figure;
plot(rad2deg(lon_true), rad2deg(lat_true), 'b', 'LineWidth', 2); hold on;
plot(rad2deg(lon_hat),  rad2deg(lat_hat),  'r--', 'LineWidth', 2);
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('Ground Track (True vs Estimated)');
legend('True', 'Estimated');
grid on;


t_sigma = t(1:length(sigma));

figure;
plot(t_sigma, rad2deg(sigma), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Bank Angle (deg)');
title('Bank Angle History');
grid on;

lat_f = lat_true(end);
lon_f = lon_true(end);

[downrange, crossrange] = compute_range_errors( ...
    r_true(end,:)', params.lat_target, params.lon_target, params.R_e);

total_error = landing_error(lat_f, lon_f, ...
    params.lat_target, params.lon_target, params.R_e);

fprintf('\n===== MPC + EKF RESULT =====\n');
fprintf('Landing Latitude (deg): %.4f\n', rad2deg(lat_f));
fprintf('Landing Longitude (deg): %.4f\n', rad2deg(lon_f));
fprintf('Downrange (m): %.2f\n', downrange);
fprintf('Crossrange (m): %.2f\n', crossrange);
fprintf('Total Error (m): %.2f\n', total_error);



fprintf('\n===== MONTE CARLO STATS =====\n');
fprintf('Mean Error      : %.2f m\n', mean(mc.total_error));
fprintf('Std Dev Error   : %.2f m\n', std(mc.total_error));
fprintf('Median Error    : %.2f m\n', median(mc.total_error));
fprintf('Max Error       : %.2f m\n', max(mc.total_error));
fprintf('Min Error       : %.2f m\n', min(mc.total_error));

fprintf('Mean Downrange  : %.2f m\n', mean(mc.downrange));
fprintf('Mean Crossrange : %.2f m\n', mean(mc.crossrange));


fprintf('\n===== 3-SIGMA STATS =====\n');
fprintf('3-sigma Downrange  : %.2f m\n', 3*std(mc.downrange));
fprintf('3-sigma Crossrange : %.2f m\n', 3*std(mc.crossrange));
fprintf('3-sigma Total Err  : %.2f m\n', 3*std(mc.total_error));



p_q    = mean(mc.q_peak > params.q_max);
p_g    = mean(mc.g_peak > params.g_max);
p_heat = mean(mc.heat_peak > params.heat_max);

fprintf('\n===== CONSTRAINT VIOLATION PROBABILITY =====\n');
fprintf('P(q > q_max)       : %.3f\n', p_q);
fprintf('P(g > g_max)       : %.3f\n', p_g);
fprintf('P(heat > heat_max) : %.3f\n', p_heat);



figure;
histogram(mc.total_error/1000, 15);
xlabel('Landing Error (km)');
ylabel('Frequency');
title('Monte Carlo Landing Error Distribution');
grid on;


figure;
scatter(mc.downrange/1000, mc.crossrange/1000, 30, 'filled');
xlabel('Downrange Error (km)');
ylabel('Crossrange Error (km)');
title('Landing Dispersion');
grid on;
axis equal;


figure;
scatter(mc.downrange/1000, mc.crossrange/1000, 30, 'filled'); hold on;
xlabel('Downrange Error (km)');
ylabel('Crossrange Error (km)');
title('Landing Dispersion with 3\sigma Ellipse');
grid on;
axis equal;

mu = [mean(mc.downrange); mean(mc.crossrange)];
Sigma = cov(mc.downrange, mc.crossrange);

theta = linspace(0, 2*pi, 200);
circle = [cos(theta); sin(theta)];

[U,S,~] = svd(Sigma);
ellipse = U * sqrt(S) * circle * 3;

plot((ellipse(1,:) + mu(1))/1000, (ellipse(2,:) + mu(2))/1000, ...
     'r', 'LineWidth', 2);
legend('Samples','3\sigma ellipse');


figure;
[f, x] = ecdf(mc.total_error);
plot(x/1000, f, 'LineWidth', 2);
xlabel('Landing Error (km)');
ylabel('Cumulative Probability');
title('CDF of Landing Error');
grid on;