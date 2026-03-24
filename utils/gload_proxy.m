function ng = gload_proxy(x, sigma, params)

    r = x(1:3);
    v = x(4:6);

    if length(x) >= 11
        beta_scale = x(7);
        rho_bias   = x(8);
        w          = x(9:11);
    else
        beta_scale = 1;
        rho_bias   = 0;
        w          = zeros(3,1);
    end

    h = norm(r) - params.R_e;
    v_rel = v - w;
    V = norm(v_rel);

    if V < 1e-6
        ng = 0;
        return;
    end

    rho = params.rho0 * exp(-h / params.H);
    rho = rho * max(0.05, 1 + rho_bias);

    q = 0.5 * rho * V^2;

    v_hat = v_rel / V;

    D = - (q / (params.beta_nom * beta_scale)) * v_hat;

    lift_dir = compute_lift_direction(r, v_rel, sigma);
    L = (q / (params.beta_nom * beta_scale)) * params.LD * lift_dir;

    a_aero = D + L;
    ng = norm(a_aero) / params.g0;
end