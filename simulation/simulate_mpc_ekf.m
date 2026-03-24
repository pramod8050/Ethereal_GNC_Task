function result = simulate_mpc_ekf(x0_true, x0_hat, P0, params)
% SIMULATE_MPC_EKF
% Closed-loop simulation using:
%   - true nonlinear dynamics
%   - EKF state estimation
%   - MPC guidance using estimated state
%
% Required fields in params:
%   params.tf
%   params.R_e
%   params.est.dt
%   params.est.Q
%   params.est.R
%   params.mpc.N
%   params.mpc.dt
%   params.mpc.u_min
%   params.mpc.u_max
%   params.mpc.options
%
% Required functions:
%   dynamics_3dof
%   ekf_predict
%   ekf_update
%   mpc_cost

    % Use one common step to avoid prediction mismatch
    dt = params.mpc.dt;
    tf = params.tf;

    % Keep estimator time step consistent with MPC
    params.est.dt = dt;

    t = 0:dt:tf;
    Nsteps = length(t);

    nx = length(x0_true);

    % Storage
    x_true_hist = zeros(nx, Nsteps);
    x_hat_hist  = zeros(nx, Nsteps);
    P_hist      = zeros(nx, nx, Nsteps);
    sigma_hist  = zeros(Nsteps-1, 1);

    % Initialize
    x_true = x0_true;
    x_hat  = x0_hat;
    P      = P0;

    x_true_hist(:,1) = x_true;
    x_hat_hist(:,1)  = x_hat;
    P_hist(:,:,1)    = P;

    % Previous MPC solution for warm start
    U_prev = zeros(params.mpc.N, 1);

    for k = 1:Nsteps-1

        % Stop if reached ground
        alt = norm(x_true(1:3)) - params.R_e;
        if alt <= 0
            t = t(1:k);
            x_true_hist = x_true_hist(:,1:k);
            x_hat_hist  = x_hat_hist(:,1:k);
            P_hist      = P_hist(:,:,1:k);
            sigma_hist  = sigma_hist(1:max(k-1,1));
            break;
        end

        % =========================================================
        % MPC CONTROL (use estimated state, NOT true state)
        % =========================================================
        [sigma, U_prev] = solve_mpc_guidance(x_hat, U_prev, params);
        sigma_hist(k) = sigma;

        % =========================================================
        % TRUE DYNAMICS
        % =========================================================
        x_true = rk4_step_true(x_true, sigma, dt, params);

        % =========================================================
        % MEASUREMENT
        % Example: direct noisy position + velocity measurement
        % =========================================================
        z = [x_true(1:3); x_true(4:6)] + ...
            chol(params.est.R, 'lower') * randn(6,1);

        % =========================================================
        % EKF PREDICT
        % =========================================================
        [x_pred, P_pred] = ekf_predict(x_hat, P, sigma, params, params.est.Q);

        % =========================================================
        % EKF UPDATE
        % =========================================================
        [x_hat, P] = ekf_update(x_pred, P_pred, z, params.est.R);

        % Store
        x_true_hist(:,k+1) = x_true;
        x_hat_hist(:,k+1)  = x_hat;
        P_hist(:,:,k+1)    = P;

        % Progress print
        if mod(k,10) == 0
            fprintf(['Step %d | Alt %.1f km | Pos err %.1f m | ', ...
                     'sigma %.2f deg\n'], ...
                     k, alt/1000, norm(x_true(1:3)-x_hat(1:3)), rad2deg(sigma));
        end
    end

    % Output
    result.t      = t(:);
    result.x_true = x_true_hist.';
    result.x_hat  = x_hat_hist.';
    result.P      = P_hist;
    result.sigma  = sigma_hist(:);
end


% ========================================================================
% MPC SOLVER
% ========================================================================
function [sigma, U_opt] = solve_mpc_guidance(x_hat, U_prev, params)
% Solve receding-horizon MPC and return first bank command

    N = params.mpc.N;

    % Warm start from shifted previous solution
    if isempty(U_prev) || length(U_prev) ~= N
        U0 = zeros(N,1);
        U_prev = U0;
    else
        U0 = [U_prev(2:end); U_prev(end)];
    end

    % Reference for consistency penalty, if used in mpc_cost
    params.mpc.U_ref = U_prev;

    % Bounds
    lb = params.mpc.u_min * ones(N,1);
    ub = params.mpc.u_max * ones(N,1);

    % Solve
    U_opt = fmincon( ...
        @(U) mpc_cost(U, x_hat, params), ...
        U0, ...
        [], [], [], [], ...
        lb, ub, ...
        [], ...
        params.mpc.options);

    % Apply first control
    sigma = U_opt(1);
end


% ========================================================================
% TRUE DYNAMICS PROPAGATION
% ========================================================================
function x_next = rk4_step_true(x, u, dt, params)
% RK4 step for truth model

    k1 = dynamics_3dof(0, x,             u, params, "nominal");
    k2 = dynamics_3dof(0, x + 0.5*dt*k1, u, params, "nominal");
    k3 = dynamics_3dof(0, x + 0.5*dt*k2, u, params, "nominal");
    k4 = dynamics_3dof(0, x + dt*k3,     u, params, "nominal");

    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end