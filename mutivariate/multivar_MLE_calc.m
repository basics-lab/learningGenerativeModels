clear;clc;
rng(0)
d = 10;
n = 50000;
A = eye(d) + randn(d,d)./(2*d);
sigma = A*A';
mu = 1 + zeros(d,1);
U_star = inv(sigma)/2;
V_star = -2*U_star*mu;
Y = max(0, mvnrnd(mu,sigma,n)');
S = Y > 0;
num_untrunc = sum(S, 2);
mu0 = zeros(d,1);
Ycent  = Y - mu0;
Ycent(~S) = 0;
var = sum(Ycent.^2,2)./num_untrunc;
U0 = diag(var.^(-1))/2;
V0 = -2*U0*mu0;
[U_hat,  V_hat] = projSGD_multivar(Y,U_star,V_star);
tik

toc

norm(U_hat - U_star, 'fro')
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