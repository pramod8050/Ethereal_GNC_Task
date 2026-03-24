function [lat_f, lon_f, x_hist] = simulate_constant_bank(x0, sigma, params)
    dt = params.mpc.dt;
    t = 0:dt:params.tf;

    x = x0;
    nx = length(x0);
    x_hist = zeros(length(t), nx);
    x_hist(1,:) = x';

    k_end = 1;
    for k = 1:length(t)-1
        alt = norm(x(1:3)) - params.R_e;
        if alt <= 0
            k_end = k;
            break;
        end

        k1 = dynamics_3dof(0, x,             sigma, params, "nominal");
        k2 = dynamics_3dof(0, x + 0.5*dt*k1, sigma, params, "nominal");
        k3 = dynamics_3dof(0, x + 0.5*dt*k2, sigma, params, "nominal");
        k4 = dynamics_3dof(0, x + dt*k3,     sigma, params, "nominal");

        x = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        x_hist(k+1,:) = x';
        k_end = k+1;
    end

    x_hist = x_hist(1:k_end,:);
    [lat_arr, lon_arr] = ecef_to_latlon(x_hist(:,1:3));
    lat_f = lat_arr(end);
    lon_f = lon_arr(end);
end