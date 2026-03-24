function mc_results = monte_carlo_ekf(x0, params, Nsim)

nx = length(x0);

% Storage
downrange_all   = zeros(Nsim,1);
crossrange_all  = zeros(Nsim,1);
total_error_all = zeros(Nsim,1);

q_peak_all      = zeros(Nsim,1);
g_peak_all      = zeros(Nsim,1);
heat_peak_all   = zeros(Nsim,1);

fprintf('\n===== MONTE CARLO START (%d runs) =====\n', Nsim);

rng('shuffle');

for i = 1:Nsim

    fprintf('Run %d / %d\n', i, Nsim);

    %% === TRUE INITIAL STATE ===
    x0_true = x0;

    % Position uncertainty (±1 km)
    x0_true(1:3) = x0_true(1:3) + 1000*randn(3,1);

    % Velocity uncertainty (±10 m/s)
    x0_true(4:6) = x0_true(4:6) + 10*randn(3,1);

    % Density bias
    x0_true(8) = 0.1 * randn;

    % Wind (m/s)
    x0_true(9:11) = 50 * randn(3,1);

    %% === INITIAL ESTIMATE ===
    x0_hat = x0;

    x0_hat(1:3) = x0_hat(1:3) + 1000*randn(3,1);
    x0_hat(4:6) = x0_hat(4:6) + 10*randn(3,1);

    %% === INITIAL COVARIANCE ===
    P0 = diag([ ...
        1e6 * ones(3,1); ...
        1e2 * ones(3,1); ...
        0.1; ...
        0.1; ...
        100 * ones(3,1) ...
    ]);

    %% === RESET MPC MEMORY ===
    clear guidance_mpc

    %% === RUN SIMULATION ===
    result = simulate_mpc_ekf(x0_true, x0_hat, P0, params);

    x_hist = result.x_true;
    sigma_hist = result.sigma;

    %% === FINAL POSITION ===
    r_final = x_hist(end,1:3)';

    [lat_f, lon_f] = ecef_to_latlon(r_final');  % ensure row vector

    %% === LANDING ERRORS ===
    [downrange, crossrange] = compute_range_errors( ...
        r_final, params.lat_target, params.lon_target, params.R_e);

    total_error = landing_error(lat_f, lon_f, ...
        params.lat_target, params.lon_target, params.R_e);

    downrange_all(i)   = downrange;
    crossrange_all(i)  = crossrange;
    total_error_all(i) = total_error;

    %% === PEAK CONSTRAINT TRACKING ===
    q_peak = 0;
    g_peak = 0;
    heat_peak = 0;

    for k = 1:size(x_hist,1)

        xk = x_hist(k,:)';

        qk = dynamic_pressure(xk, params);
        heatk = heat_rate_proxy(xk, params);

        if k <= length(sigma_hist)
            sigmak = sigma_hist(k);
        else
            sigmak = sigma_hist(end);
        end

        gk = gload_proxy(xk, sigmak, params);

        q_peak    = max(q_peak, qk);
        g_peak    = max(g_peak, gk);
        heat_peak = max(heat_peak, heatk);
    end

    q_peak_all(i)    = q_peak;
    g_peak_all(i)    = g_peak;
    heat_peak_all(i) = heat_peak;

end

%% === STORE RESULTS ===
mc_results.downrange   = downrange_all;
mc_results.crossrange  = crossrange_all;
mc_results.total_error = total_error_all;

mc_results.q_peak    = q_peak_all;
mc_results.g_peak    = g_peak_all;
mc_results.heat_peak = heat_peak_all;

fprintf('===== MONTE CARLO COMPLETE =====\n');

end