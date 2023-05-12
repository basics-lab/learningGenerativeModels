function plot_mle_vs_n(dist_idx, n_dist_mc, n_mc, seed)
%% Setup
distributions = ["Normal", "Exponential", "Cauchy"];
distribution = distributions(dist_idx);
d = 10;
rng(seed)
N=10; % number of points
n_space = ceil(logspace(2, 5, N)); %number of samples
w_star = randn(d,1);
sigma = 1;
[u_star,v_star] = parameter_tf(w_star, sigma);
%%
dist_v_n = zeros(N,3);
for j = 1:N
    n = n_space(j);
    dist = zeros(n_mc, 1);
    fprintf("Working on n=%i\n",n)
    tic;
    parfor i = 1:n_mc
        % Get new samples
        if distribution == "Cauchy"
            X = trnd(1, d, n);
        elseif distribution == "Normal"
            X = randn(d, n);
        elseif distribution == "Exponential"
            X = (2*(rand(d,n) > 0.5) - 1).*exprnd(1,d,n);
        end
        y = (w_star'*X)' + sigma*randn(n,1);
        u0 = randn(d,1);
        v0 = 0.8;
        tic;
        %[u_hat_cvx, v_hat_cvx] = cvx_solve(X,y);
        [u_hat, v_hat, res, iter] = projGD_MLE(y, X, u0, v0);
        fprintf("MLE Compute time: %.3f, in %d iterations, res=%d\n", toc, iter, res);
        tic;
        if distribution == "Cauchy"
            X_sample = trnd(1, d, n_dist_mc);
        elseif distribution == "Normal"
            X_sample = randn(d, n_dist_mc);
        elseif distribution == "Exponential"
            X_sample = (2*(rand(d,n_dist_mc) > 0.5) - 1).*exprnd(1,d,n_dist_mc);
        end
        dist(i) = distance_metrics(u_hat, v_hat, u_star, v_star, X_sample);
        fprintf('distance compute time=%.3f\n', toc);
    end
    time = toc;
    fprintf("Total time = %.2fsec\n",time)
    dist_v_n(j,:) = mean(dist);
end
filename = sprintf('%s_seed%d_%s.mat', distribution, seed, datestr(now,'HH_MM_SS_FFF'));
save(filename);
%% Helper Functions
function [u,v] = parameter_tf(w,s)
    u=w/s; v=1/s;
end
end