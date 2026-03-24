function A = numerical_jacobian(fun, x)
% NUMERICAL_JACOBIAN
% Finite-difference Jacobian of vector function fun(x)

    fx = fun(x);
    n = length(x);
    m = length(fx);

    A = zeros(m, n);
    eps_fd = 1e-6;

    for i = 1:n
        dx = zeros(n,1);
        dx(i) = eps_fd;

        A(:,i) = (fun(x + dx) - fx) / eps_fd;
    end
end