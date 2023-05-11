clear;clc;
rng(0)
d = 10;
n = 5000;
A = eye(d) + randn(d,d)./(2*d);
sigma = A*A';
mu = zeros(d,1);
Y = max(0, mvnrnd(mu,sigma,n)');
U0 = inv(sigma)/2;
V0 = -2*U0*mu;
[U_hat,  V_hat] = projSGD_multivar(Y,U0,V0);
% for i=1:5
%     [gradu, gradv] = grad_lik_function(Y,U,V)
% end
% for i=1:5
%     lik = lik_func_multivar(Y,U,V)
% end
% %% Compute the gradient via finite differences
% eps = 1e-2;
% dV = zeros(d,1);
% dU = zeros(d,d);
% for i=1:d
%     dV(i) = eps;
%     emp_gradv(i) = (lik_func_multivar(Y,U,V+dV) - lik_func_multivar(Y,U,V-dV))/(2*eps);
%     dV(i) = 0;
% end