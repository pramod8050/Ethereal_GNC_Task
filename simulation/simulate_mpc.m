function result = simulate_mpc(x0, params)

dt = params.mpc.dt;
tf = params.tf;

t = 0:dt:tf;
Nsteps = length(t);

nx = length(x0);

x = zeros(nx, Nsteps);
sigma_hist = zeros(1, Nsteps-1);

x(:,1) = x0;

for k = 1:Nsteps-1

    xk = x(:,k);

    % Stop if reached ground
    alt = norm(xk(1:3)) - params.R_e;
    if alt <= 0
        t = t(1:k);
        x = x(:,1:k);
        sigma_hist = sigma_hist(1:max(k-1,1));
        break;
    end

    if mod(k,5) == 0
        fprintf('Step %d / %d | Time: %.1f s | Altitude: %.1f km\n', ...
            k, Nsteps, t(k), alt/1000);
    end

    % MPC control
    sigma = guidance_mpc(t(k), xk, params);
    sigma_hist(k) = sigma;

    % RK4 propagation
    k1 = dynamics_3dof(0, xk, sigma, params, "nominal");
    k2 = dynamics_3dof(0, xk + 0.5*dt*k1, sigma, params, "nominal");
    k3 = dynamics_3dof(0, xk + 0.5*dt*k2, sigma, params, "nominal");
    k4 = dynamics_3dof(0, xk + dt*k3, sigma, params, "nominal");

    x(:,k+1) = xk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

result.t = t(:);
result.x = x.';
result.sigma_hist = sigma_hist(:);
end