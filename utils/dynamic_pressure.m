function q = dynamic_pressure(x, params)

    r = x(1:3);
    v = x(4:6);

    if length(x) >= 11
        rho_bias = x(8);
        w = x(9:11);
    else
        rho_bias = 0;
        w = zeros(3,1);
    end

    h = norm(r) - params.R_e;
    v_rel = v - w;
    V = norm(v_rel);

    rho = params.rho0 * exp(-h / params.H);
    rho = rho * max(0.05, 1 + rho_bias);

    q = 0.5 * rho * V^2;
end