function lift_dir = compute_lift_direction(r, v_rel, sigma)

v_hat = v_rel / norm(v_rel);
r_hat = r / norm(r);

% FIXED: correct orientation
h = cross(r_hat, v_hat);
h_hat = h / norm(h);

% Nominal lift direction
n0 = cross(h_hat, v_hat);

% Bank rotation
lift_dir = cos(sigma) * n0 + sin(sigma) * h_hat;

% Normalize
lift_dir = lift_dir / norm(lift_dir);

end