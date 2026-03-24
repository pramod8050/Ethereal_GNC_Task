function [x_pred, P_pred] = ekf_predict(x_hat, P, u, params, Q)
% EKF prediction step using one-step RK4 propagation and numerical Jacobian

    dt = params.est.dt;

    % Nonlinear state propagation
    f_disc = @(x) rk4_step(x, u, dt, params);

    x_pred = f_disc(x_hat);

    % Numerical Jacobian of discrete dynamics
    F = numerical_jacobian(f_disc, x_hat);

    % Covariance prediction
    P_pred = F * P * F' + Q;
end


function x_next = rk4_step(x, u, dt, params)
    k1 = dynamics_3dof(0, x,             u, params, "nominal");
    k2 = dynamics_3dof(0, x + 0.5*dt*k1, u, params, "nominal");
    k3 = dynamics_3dof(0, x + 0.5*dt*k2, u, params, "nominal");
    k4 = dynamics_3dof(0, x + dt*k3,     u, params, "nominal");

    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end