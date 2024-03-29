function plot_mle_vs_d(dist_idx, n_dist_mc, n_mc, n, seed, mean, var)
%% Setup
if getenv('USER') == ""
    parpool(str2num(getenv('SLURM_CPUS_ON_NODE')));
    fprintf("ON SLURM, CREATING A BIGGER POOL\n")
end
distributions = ["Normal", "Exponential", "Cauchy", "Mixture"];
distribution = distributions(dist_idx);
rng(seed)
N=10; % number of points
d_space = ceil(logspace(log10(5),log10(50),N)); %number of samples
sigma = 1;
filename = sprintf('%s_seed%d_n%i_%s.mat', distribution, seed,n ,datestr(now,'HH_MM_SS_FFF'));
%%
iter_v_n = zeros(N,n_mc);
res_v_n = zeros(N,n_mc);
dist_v_n = zeros(N,n_mc);
p=0.01;
for j = 1:N
    d = d_space(j);
    w_star = ones(d,1);
    [u_star,v_star] = parameter_tf(w_star, sigma);
    dist = zeros(n_mc, 1);
    res = zeros(n_mc,1);
    iter = zeros(n_mc,1);
    fprintf("Working on d=%i\n",d)
    tic;
    parfor i = 1:n_mc
        % Get new samples
        if distribution == "Cauchy"
            X = trnd(1, d, n);
        elseif distribution == "Normal"
            X = sqrt(var)*randn(d, n) + mean;
        elseif distribution == "Exponential"
            X = (2*(rand(d,n) > 0.5) - 1).*exprnd(1,d,n);
        elseif distribution == "Mixture"
            X = randn(d, n) + mean*(rand(d,n) < p);
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
        elseif distribution == "Mixture"
            X_sample = randn(d, n) + mean*(rand(d,n) < p);
        end
        dist(i) = distance_metrics(u_hat, v_hat, u_star, v_star, X_sample);
        
        fprintf('distance compute time=%.3f\n', toc);
    end
    time = toc;
    fprintf("Total time = %.2fsec\n",time)
    dist_v_n(j,:) = dist;
    res_v_n(j,:) = res;
    iter_v_n(j,:) = iter;
    save(filename);
end
%% Helper Functions
function [u,v] = parameter_tf(w,s)
    u=w/s; v=1/s;
end
end
