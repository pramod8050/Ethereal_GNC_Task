function [x_upd, P_upd, K, innov] = ekf_update(x_pred, P_pred, z, R)
% EKF measurement update for z = [r; v]

    nx = length(x_pred);

    % Measurement model: h(x) = [r; v]
    h = [x_pred(1:3); x_pred(4:6)];

    % Measurement Jacobian
    H = zeros(6, nx);
    H(1:3,1:3) = eye(3);
    H(4:6,4:6) = eye(3);

    innov = z - h;

    S = H * P_pred * H' + R;
    K = P_pred * H' / S;

    x_upd = x_pred + K * innov;

    I = eye(nx);
    P_upd = (I - K*H) * P_pred;
end