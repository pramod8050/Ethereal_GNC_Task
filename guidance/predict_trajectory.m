function X = predict_trajectory(x0, U, params)

    N  = length(U);
    dt = params.mpc.dt;
    nx = length(x0);

    X = zeros(nx, N+1);
    X(:,1) = x0;

    xk = x0;

    for k = 1:N
        uk = U(k);

        k1 = dynamics_3dof(0, xk,              uk, params, "nominal");
        k2 = dynamics_3dof(0, xk + 0.5*dt*k1, uk, params, "nominal");
        k3 = dynamics_3dof(0, xk + 0.5*dt*k2, uk, params, "nominal");
        k4 = dynamics_3dof(0, xk + dt*k3,     uk, params, "nominal");

        xk = xk + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

        X(:,k+1) = xk;
    end
end