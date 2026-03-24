function result = simulate_single(x0, params, mode)

    tspan = [0 params.tf];

    u_func = @(t, x) guidance_mpc(t, x, params);

    options = odeset('Events', @(t,x) event_ground(t,x,params));

    [t, x] = ode45(@(t,x) dynamics_3dof(t, x, u_func(t,x), params, mode), ...
                   tspan, x0, options);

    result.t = t;
    result.x = x;
end