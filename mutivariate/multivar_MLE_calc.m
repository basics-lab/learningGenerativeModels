clear;clc;
rng(1)
d = 3;
n = 10000;
kappa = 3; %sqrt condition number
A = eye(d) + randn(d,d)./(2*d);
sigma = A*A';
[L,D] = eig(sigma);
D = diag(min(max(1/kappa, diag(D)), kappa));
D(1,1) = kappa;
D(2,2) = 1/kappa;
sigma = L*D*L';
mu = 0.5 + zeros(d,1);
U_star = inv(sigma)/2;
V_star = -2*U_star*mu;
Y = max(0, mvnrnd(mu,sigma,n)');
tic;
[sigma_hat, b_hat] = main(Y);
%sigma_hat = sigma;
%b_hat = mu;
U_hat_wu = inv(sigma_hat)/2;
V_hat_wu = -2*U_hat_wu*b_hat;
time = toc;
V0 = -ones(d,1);
U0 = eye(d)/2;
fprintf("Shanshan: runtime=%.3f\n", time);
[U_hat,  V_hat, U_final, V_final] = projSGD_multivar(Y,U0,V0);
%% Metrics
norm(U_hat - U_star, 'fro')
norm(inv(sigma_hat)/2 - U_star, 'fro')
n_sample = 4000;
Y = max(0, mvnrnd(mu,sigma,n_sample)');
% KL computation
fprintf("Computing KL\n")
KL = zeros(1,n_sample);
KL_wu = zeros(1,n_sample);
for i=1:n_sample
    KL(i) = -multivar_LLR(Y(:,i), U_hat, V_hat, U_star, V_star);
    KL_wu(i) = -multivar_LLR(Y(:,i), U_hat_wu, V_hat_wu, U_star, V_star);
end
d1 = sqrt(0.5*sum(KL)/n_sample);
d2 = sqrt(0.5*sum(KL_wu)/n_sample);
fprintf("KL of ours: %.3f\nKL of Wu: %.3f\n", d1,d2);
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