function sigma = guidance_mpc(t, x, params)
% GUIDANCE_MPC
% Returns first control action from finite-horizon MPC

    persistent U_prev sigma_prev

    N = params.mpc.N;

    % Initial guess
    if isempty(U_prev) || length(U_prev) ~= N
        U0 = deg2rad(20) * ones(N,1);
        U_prev = U0;
    else
        U0 = [U_prev(2:end); U_prev(end)];
    end

    if isempty(sigma_prev)
        sigma_prev = U0(1);
    end

    % Bounds
    lb = params.mpc.u_min * ones(N,1);
    ub = params.mpc.u_max * ones(N,1);

    % Pass warm-start reference into cost
    params.mpc.U_ref = U_prev;

    % Pass previous applied control into nonlinear constraint
    params.mpc.sigma_prev = sigma_prev;

    % Cost handle
    J = @(U) mpc_cost(U, x, params);

    % Nonlinear constraints
    nonlcon = @(U) mpc_nonlcon(U, params);

    % Solve
    try
        U_opt = fmincon(J, U0, [], [], [], [], lb, ub, nonlcon, params.mpc.options);
    catch
        U_opt = U0;
    end

    % First move
    sigma = U_opt(1);

    % Final clamp
    sigma = max(min(sigma, params.mpc.u_max), params.mpc.u_min);

    % Save for next solve
    U_prev = U_opt;
    sigma_prev = sigma;
end

function [c, ceq] = mpc_nonlcon(U, params)
% Nonlinear inequality constraints for bank-rate limiting
% c(U) <= 0

    ceq = [];
    c = [];

    if isfield(params, 'sigma_dot_max') && ~isempty(params.sigma_dot_max)
        rate_max = params.sigma_dot_max;
        dt = params.mpc.dt;

        % First move relative to previously applied control
        if isfield(params.mpc, 'sigma_prev')
            c0 = abs(U(1) - params.mpc.sigma_prev)/dt - rate_max;
            c = [c; c0];
        end

        % Remaining moves inside horizon
        if length(U) > 1
            dU = diff(U)/dt;
            c_rate = abs(dU) - rate_max;
            c = [c; c_rate(:)];
        end
    end
end