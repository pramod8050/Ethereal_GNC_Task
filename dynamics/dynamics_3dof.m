function dx = dynamics_3dof(t, x, u, params, mode)

% States
r = x(1:3);
v = x(4:6);

% Optional augmented states
beta_scale = x(7);
rho_bias   = x(8);
w          = x(9:11);

% Relative velocity
v_rel = v - w;
v_rel_norm = norm(v_rel);

% Altitude
h = norm(r) - params.R_e;

% Density
rho = params.rho0 * exp(-h / params.H) * (1 + rho_bias);

% Dynamic pressure
q = 0.5 * rho * v_rel_norm^2;

% Drag
D = - (q / (params.beta_nom * beta_scale)) * (v_rel / v_rel_norm);

% Lift direction (simplified)
% You will refine this later
lift_dir = compute_lift_direction(r, v_rel, u);

L = (q / (params.beta_nom * beta_scale)) * params.LD * lift_dir;

aero = D + L;

% Gravity
g = -params.mu * r / norm(r)^3;

% Rotation effects
omega = params.omega;
coriolis = -2 * cross(omega, v);
centrifugal = -cross(omega, cross(omega, r));

% Dynamics
dr = v;
dv = g + aero + coriolis + centrifugal;

% Augmented states (random walk)
dbeta = 0;
drho  = 0;
dw    = zeros(3,1);

dx = [dr; dv; dbeta; drho; dw];

end