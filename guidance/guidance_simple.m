function sigma = guidance_simple(t, x, params)

r_current = x(1:3);

[downrange, crossrange] = compute_range_errors( ...
    r_current, params.lat_target, params.lon_target, params.R_e);
k = params.k_guidance;


sigma = -k * crossrange;

% saturation
sigma = max(min(sigma, params.sigma_max), -params.sigma_max);



end