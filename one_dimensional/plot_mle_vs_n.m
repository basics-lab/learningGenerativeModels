function plot_mle_vs_n(dist_idx, n_dist_mc, n_mc, d, seed, mean, var)
%% Setup
if getenv('USER') == "justinkang"
    parpool(str2num(getenv('SLURM_CPUS_ON_NODE')));
    fprintf("ON SLURM, CREATING A BIGGER POOL\n")
end
distributions = ["Normal", "Exponential", "Cauchy"];
distribution = distributions(dist_idx);
rng(seed)
N=10; % number of points
n_space = ceil(logspace(2, 5, N)); %number of samples
w_star = ones(d,1);
sigma = 1;
[u_star,v_star] = parameter_tf(w_star, sigma);
%%
iter_v_n = zeros(N,n_mc);
res_v_n = zeros(N,n_mc);
dist_v_n = zeros(N,n_mc);
for j = 1:N
    n = n_space(j);
    dist = zeros(n_mc, 1);
    res = zeros(n_mc,1);
    iter = zeros(n_mc,1);
    fprintf("Working on n=%i\n",n)
    tic;
    parfor i = 1:n_mc
        % Get new samples
        if distribution == "Cauchy"
            X = trnd(1, d, n);
        elseif distribution == "Normal"
            X = sqrt(var)*randn(d, n) + mean;
        elseif distribution == "Exponential"
            X = (2*(rand(d,n) > 0.5) - 1).*exprnd(1,d,n);
        end
        y = (w_star'*X)' + sigma*randn(n,1);
        u0 = randn(d,1);
        v0 = 0.8;
        tic;
        %[u_hat_cvx, v_hat_cvx] = cvx_solve(X,y);
        [u_hat, v_hat, res(i), iter(i)] = projGD_MLE(y, X, u0, v0);
        fprintf("MLE Compute time: %.3f, in %d iterations, res=%d\n", toc, iter(i), res(i));
        tic;
        if distribution == "Cauchy"
            X_sample = trnd(1, d, n_dist_mc);
        elseif distribution == "Normal"
            X_sample = sqrt(var)*randn(d, n_dist_mc) + mean;
        elseif distribution == "Exponential"
            X_sample = (2*(rand(d,n_dist_mc) > 0.5) - 1).*exprnd(1,d,n_dist_mc);
        end
        dist(i) = distance_metrics(u_hat, v_hat, u_star, v_star, X_sample);
        
        fprintf('distance compute time=%.3f\n', toc);
    end
    time = toc;
    fprintf("Total time = %.2fsec\n",time)
    dist_v_n(j,:) = dist;
    res_v_n(j,:) = res;
    iter_v_n(j,:) = iter;
end
filename = sprintf('%s_seed%d_d%i_%s.mat', distribution, seed,d ,datestr(now,'HH_MM_SS_FFF'));
save(filename);
%% Helper Functions
function [u,v] = parameter_tf(w,s)
    u=w/s; v=1/s;
end
end