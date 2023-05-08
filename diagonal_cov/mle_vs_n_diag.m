clear;clc;
rng(0);
d = 3;
n= 10000;
b = 5*ones(d,1);
sigma = 0.5 + exprnd(0.5, d,1);
u_star = b./sigma;
v_star = 1./sigma;
u0 = b + randn(d,1);
v0 = 0.5 + exprnd(0.5, d,1);
Y = max(0, ((sigma*ones(1,n)).*randn(d,n)) + b);
%% Blackbox Optimizer (Neldar-Mead)
l = @(in) -lik_func_diag(Y, in(1:d), in(d+1:2*d)); %Objective
tic
uv_nm = fminsearch(l, [u_star; v_star]);
t_nm = toc;
%% Quasi-Newton
options = optimoptions('fminunc','SpecifyObjectiveGradient',true, "Algorithm","quasi-newton");
l = @(x) fminobj(x, Y); %Objective
tic
uv_qn = fminunc(l,[u0; v0],options);
t_qn = toc;
%% Trust-Region
options = optimoptions('fminunc','SpecifyObjectiveGradient',true, "Algorithm","trust-region");
l = @(x) fminobj(x, Y); %Objective
tic
uv_tr = fminunc(l,[u0; v0],options);
t_tr = toc;
%% PGD with line Search
tic
[u_hat, v_hat] = projGD_MLE_diag(Y, u0, v0);
t_pgd = toc;
%% Shanshan 
[sigma_hat_wu, b_hat_wu] = main_diag(Y, u0./v0,1./v0);
sigma_hat_wu = diag(sigma_hat_wu).^(0.5);
%% Negative Log Likelihood
fprintf("NLL of Truth: %.3f\n", l([u_star; v_star]))
fprintf("NLL of NM Solution %.3f in %.4f sec\n", l(uv_nm), t_nm)
fprintf("NLL of QN Solution %.3f in %.4f sec\n", l(uv_qn), t_qn)
fprintf("NLL of TR Solution %.3f in %.4f sec\n", l(uv_tr), t_tr)
fprintf("NLL of PGD Solution %.3f in %.4f sec\n", l([u_hat; v_hat]), t_pgd)
%% Distance Metrics
X_sample = randn(d, n);
dist = dist_prod_distr(u_hat,v_hat,u_star,v_star);
fprintf("TV of MLE: %.3f\n", dist(1));
fprintf("KL of MLE: %.3f\n", dist(2));
fprintf("TV Upper bound of MLE: %.3f\n", dist(3));

function [f, GRAD] = fminobj(x, Y)
    [d,~] = size(Y);
    f = -lik_func_diag(Y, x(1:d), x(d+1:2*d));
    GRAD = grad_lik_func_diag(Y, x(1:d), x(d+1:2*d));
end

function dist = dist_prod_distr(u1,v1,u2,v2)
    dist(1) = compute_TV(v1.^(-1), u1./v1, v2.^(-1), u2./v2);
    dist(2) = (0.5*KL_div(u1./v1, diag(v1.^(-2)), u2./v2, diag(v2.^(-2)))).^(0.5);
    dist(3) = tv_tringle(v1.^(-1), u1./v1, v2.^(-1), u2./v2);
end

function dist = tv_tringle(s1,b1,s2,b2)
    dist = sum(abs(normcdf(-b1./s1) - normcdf(-b2./s2)));
    diff = @(x) abs(normpdf(x, b1, s1) - normpdf(x, b2, s2));
    dist  = dist + sum(integral(diff, 0, max([b1;b2]) + 6*max([s1;s2]), 'ArrayValued', true));
end