function params = config()

params.mu = 3.986e14;             % Earth gravitational parameter
params.R_e = 6371e3;              % Earth radius
params.omega = [0; 0; 7.292e-5];  % Earth rotation

% Aerodynamics
params.beta_nom = 1000;           
params.LD = 0.3;

% Atmosphere
params.rho0 = 1.225;
params.H = 7500;

% Gravity constant
params.g0 = 9.81;

% Simulation
params.dt = 0.1;
params.tf = 1000;

% Target
params.lat_target = deg2rad(0.3);
params.lon_target = deg2rad(9.2);

% Guidance
params.k_guidance = 1e-5;

% Heat model coefficient
params.k_heat = 1e-7;

% Path limits
params.q_max = 8e4;                 
params.g_max = 8;                   
params.heat_max = 2e6;

% Bank limits
params.sigma_max = deg2rad(70);
params.sigma_dot_max = deg2rad(5);

%% =============================
% MPC SETTINGS
% =============================
params.mpc.N  = 40;
params.mpc.dt = 2.0;

% Control limits
params.mpc.u_min = -params.sigma_max;
params.mpc.u_max =  params.sigma_max;

% Terminal cost
params.mpc.w_terminal_dr = 50;
params.mpc.w_terminal_cr = 200;

% Stage cost
params.mpc.w_stage = 0.0;
params.mpc.w_dr    = 300;
params.mpc.w_cr    = 150;

% Control regularization
params.mpc.w_u  = 1e-6;
params.mpc.w_du = 1e-2;

% MPC consistency
params.mpc.w_track = 1e-1;

% Heat / soft constraints
params.mpc.w_heat = 0.0;
params.mpc.use_soft_constraints = false;
params.mpc.w_soft = 1e-6;

% Optimizer
params.mpc.options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'none', ...
    'MaxIterations', 50, ...
    'MaxFunctionEvaluations', 5000, ...
    'OptimalityTolerance', 1e-3, ...
    'StepTolerance', 1e-3);

%% =============================
% EKF SETTINGS
% =============================
params.est.dt = params.mpc.dt;

q_r    = 1e-3;
q_v    = 1e-2;
q_beta = 1e-6;
q_rho  = 1e-6;
q_w    = 1e-3;

params.est.Q = diag([ ...
    q_r, q_r, q_r, ...
    q_v, q_v, q_v, ...
    q_beta, ...
    q_rho, ...
    q_w, q_w, q_w]);

r_pos = 10;
r_vel = 1;

params.est.R = diag([ ...
    r_pos^2, r_pos^2, r_pos^2, ...
    r_vel^2, r_vel^2, r_vel^2]);

params.mpc.w_flip = 1e3;
params.mpc.w_hN   = 1e-2;
params.mpc.w_vN   = 1e-4;

end