clear;clc;
rng(4);
d = 3;
n= 5000;
b = ones(d,1);
sigma = [5;0.2;1];
u_star = b./sigma;
v_star = 1./sigma;
u0 = b + randn(d,1);
v0 = 0.5 + exprnd(0.5, d,1);
Y = max(0, ((sigma*ones(1,n)).*randn(d,n)) + b);
l = @(in) -lik_func_diag(Y, in(1:d), in(d+1:2*d)); %Objective
%% PGD with line Search
tic
[u_hat, v_hat] = projGD_MLE_diag(Y, u0, v0);
t_pgd = toc;
%% Shanshan 
[sigma_hat_wu, b_hat_wu] = main_diag(Y, u0./v0,1./v0);
v_hat_wu = diag(sigma_hat_wu).^(-0.5);
u_hat_wu = b_hat_wu'.*v_hat_wu;
%% Negative Log Likelihood
fprintf("NLL of Truth: %.3f\n", l([u_star; v_star]))
fprintf("NLL of PGD Solution %.3f in %.4f sec\n", l([u_hat; v_hat]), t_pgd)
%% Distance Metrics
X_sample = randn(d, n);
dist_mle = dist_prod(u_hat,v_hat,u_star,v_star);
fprintf("MLE Dist from True:\n")
print_dist(dist_mle)
dist_wu = dist_prod(u_hat_wu,v_hat_wu,u_star,v_star);
fprintf("\nWu Dist from True:\n")
print_dist(dist_wu)

function [f, GRAD] = fminobj(x, Y)
    [d,~] = size(Y);
    f = -lik_func_diag(Y, x(1:d), x(d+1:2*d));
    GRAD = grad_lik_func_diag(Y, x(1:d), x(d+1:2*d));
end

function print_dist(dist)
    fprintf("TV: %.3f\n", dist(1));
    fprintf("Pinsker KL (untrunc): %.3f\n", dist(2));
    fprintf("Pinsker Reverse KL (untrunc): %.3f\n", dist(3));
    fprintf("Pinsker True KL: %.3f\n", dist(4));
    fprintf("Pinsker Reverse True KL bound of MLE: %.3f\n", dist(5));
    fprintf("TV triangle inequality MLE: %.3f\n", dist(6));
end



