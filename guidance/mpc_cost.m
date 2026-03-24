function J = mpc_cost(U, x0, params)
% MPC_COST

    % Predict trajectory
    X = predict_trajectory(x0, U, params);

    % Terminal state
    xN = X(:,end);
    rN = xN(1:3);
    vN = xN(4:6);

    % Terminal range errors
    [downrange, crossrange] = compute_range_errors( ...
        rN, params.lat_target, params.lon_target, params.R_e);

    % Terminal miss cost
    J_terminal = params.mpc.w_terminal_dr * downrange^2 + ...
                 params.mpc.w_terminal_cr * crossrange^2;

    % Terminal altitude / speed shaping
    hN = norm(rN) - params.R_e;
    VN = norm(vN);

    J_terminal_shape = 0;
    if isfield(params.mpc, 'w_hN')
        J_terminal_shape = J_terminal_shape + params.mpc.w_hN * hN^2;
    end
    if isfield(params.mpc, 'w_vN')
        J_terminal_shape = J_terminal_shape + params.mpc.w_vN * VN^2;
    end

    % Control effort
    J_u = params.mpc.w_u * sum(U.^2);

    % Smoothness
    if length(U) > 1
        dU = diff(U);
        J_du = params.mpc.w_du * sum(dU.^2);
    else
        J_du = 0;
    end

    % Warm-start / consistency tracking
    J_track = 0;
    if isfield(params.mpc, 'w_track') && params.mpc.w_track > 0 && ...
       isfield(params.mpc, 'U_ref') && ~isempty(params.mpc.U_ref)
        U_ref = params.mpc.U_ref(:);
        if length(U_ref) == length(U)
            J_track = params.mpc.w_track * sum((U - U_ref).^2);
        end
    end

    % Penalize sign changes
    J_flip = 0;
    if isfield(params.mpc, 'w_flip') && params.mpc.w_flip > 0 && length(U) > 1
        for k = 2:length(U)
            if U(k)*U(k-1) < 0
                J_flip = J_flip + params.mpc.w_flip;
            end
        end
    end

    % Stage cost
    J_stage = 0;
    if isfield(params.mpc, 'w_stage') && params.mpc.w_stage > 0
        N = length(U);
        k_start = max(1, round(0.8 * N));

        for k = k_start:N
            xk = X(:,k);
            rk = xk(1:3);

            [dr_k, cr_k] = compute_range_errors( ...
                rk, params.lat_target, params.lon_target, params.R_e);

            alpha = (k - k_start + 1) / (N - k_start + 1);

            J_stage = J_stage + alpha * params.mpc.w_stage * ...
                (params.mpc.w_dr * dr_k^2 + params.mpc.w_cr * cr_k^2);
        end
    end

    % Heat cost
    J_heat = 0;
    if isfield(params.mpc, 'w_heat') && params.mpc.w_heat > 0
        for k = 1:length(U)
            xk = X(:,k);
            heat_k = heat_rate_proxy(xk, params);
            J_heat = J_heat + params.mpc.w_heat * heat_k * params.mpc.dt;
        end
    end

    % Soft constraints
    J_soft = 0;
    if isfield(params.mpc, 'use_soft_constraints') && params.mpc.use_soft_constraints
        for k = 1:length(U)
            xk = X(:,k);

            qk    = dynamic_pressure(xk, params);
            gk    = gload_proxy(xk, U(min(k,end)), params);
            heatk = heat_rate_proxy(xk, params);

            if isfield(params, 'q_max')
                J_soft = J_soft + params.mpc.w_soft * max(0, qk - params.q_max)^2;
            end
            if isfield(params, 'g_max')
                J_soft = J_soft + params.mpc.w_soft * max(0, gk - params.g_max)^2;
            end
            if isfield(params, 'heat_max')
                J_soft = J_soft + params.mpc.w_soft * max(0, heatk - params.heat_max)^2;
            end
        end
    end

    % Total cost
    J = J_terminal + J_terminal_shape + J_stage + J_u + J_du + ...
        J_track + J_flip + J_heat + J_soft;
end