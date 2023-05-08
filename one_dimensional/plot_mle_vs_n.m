%% Setup
clear;clc;
d = 3;
N=20; % number of points
n_space = ceil(logspace(2, 5, N)); %number of samples
n_mc = 10; % number of
n_dist_mc = 1000;
w_star = randn(d,1);
sigma = 1;
[u_star,v_star] = parameter_tf(w_star, sigma);
%%
dist_v_n = zeros(N,3);
for j = 1:N
    n = n_space(j);
    dist = zeros(n_mc, 3);
    fprintf("Working on n=%i\n",n)
    tic;
    parfor i = 1:n_mc
        % Get new samples
        X = randn(d, n);
        y = (w_star'*X)' + sigma*randn(n,1);
        u0 = randn(d,1);
        v0 = 0.8;
        [u_hat, v_hat] = projGD_MLE(y, X, u0, v0);
        X_sample = randn(d, n_dist_mc);
        dist(i,:) = distance_metrics(u_hat, v_hat, u_star, v_star, X_sample);
    end
    time = toc;
    fprintf("Total time = %.2fsec\n",time)
    dist_v_n(j,:) = mean(dist,1);
end
semilogx(n_space, dist_v_n)
%% Helper Functions
function [u,v] = parameter_tf(w,s)
    u=w/s; v=1/s;
end