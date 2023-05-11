clear;clc;
rng(0)
d = 10;
n = 100;
A = eye(d) + randn(d,d)./(2*d);
sigma = A*A';
mu = zeros(d,1);
Y = max(0, mvnrnd(mu,sigma,n)');
U = inv(sigma)/2;
V = -2*U*mu;
for i=1:5
    [gradu, gradv] = grad_lik_function(Y,U,V)
end
for i=1:5
    lik = lik_func_multivar(Y,U,V)
end
%% Compute the gradient via finite differences
eps = 1e-3;
dV = zeros(d,1);
dU = zeros(d,d);
for i=1:d
    dV(i) = eps;
    emp_gradv(i) = (lik_func_multivar(Y,U,V+dV) - lik_func_multivar(Y,U,V-dV))/(2*eps);
    dV(i) = 0;
end
for i=1:d
    for j=1:d
        dU(i,j) = eps;
        emp_gradU(i,j) = (lik_func_multivar(Y,U+dU,V) - lik_func_multivar(Y,U-dU,V))/(2*eps);
        dU(i,j) = 0;
    end
end