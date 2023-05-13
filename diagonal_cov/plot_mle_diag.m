function plot_mle_diag(N,n_mc,n, seed)
d = 3;
rng(seed);
b_space = linspace(-5,5,N); %grid
sigma = ones(3,1);
verbose = true;
%% 
dist_v_b = zeros(N,2);
dist_v_b_wu = zeros(N,2);
dist_v_b_wu2 = zeros(N,2);
dist = zeros(n_mc,2,N);
iter = zeros(n_mc,N);
res = zeros(n_mc,N);
dist_wu = zeros(n_mc,2,N);
zero_sample_wu = zeros(n_mc,N);
dist_wu2 = zeros(n_mc,2,N);
zero_sample_wu2 = zeros(n_mc,N);
iter_wu = zeros(n_mc,N);
res_wu = zeros(n_mc,N);
for j = 1:N
    b = b_space(j)*ones(d,1);
    u_star = b./sigma;
    v_star = 1./sigma;
    tic;
    parfor i = 1:n_mc
        % Get new samples
        Y = max(0, ((sigma*ones(1,n)).*randn(d,n)) + b);
        u0 = b(1) + randn(d,1);
        v0 = 0.5 + exprnd(0.5, d,1);
        % Our Algorithm
        [u_hat, v_hat, res(i,j), iter(i,j)] = projGD_MLE_diag(Y, u0, v0);
        dist(i,:,j) = dist_prod_distr(u_hat, v_hat, u_star, v_star);
        % Wu et al. Modified
        [sigma_hat_wu, b_hat_wu, res_wu(i,j), iter_wu(i,j), zero_sample_wu(j,i)] = main_diag(Y, u0, v0);
        v_hat_wu = diag(sigma_hat_wu).^(0.5);
        u_hat_wu = b_hat_wu'./v_hat_wu;
        dist_wu(i,:,j) = dist_prod_distr(u_hat_wu, v_hat_wu, u_star, v_star);
        % Wu et al. unmodified
        if b > 0
            [sigma_hat_wu2, b_hat_wu2, zero_sample_wu2(j,i)] = main_diag_PGD(Y);
            v_hat_wu2 = diag(sigma_hat_wu2).^(0.5);
            u_hat_wu2 = b_hat_wu2'./v_hat_wu2;
            dist_wu2(i,:,j) = dist_prod_distr(u_hat_wu2, v_hat_wu2, u_star, v_star);
        end
        if verbose
            fprintf("Opt took t=%i iterations with residual=%f\n",iter(i,j),res(i,j))
            fprintf("Wu Opt took t_avg=%.1f iterations with avg. residual=%f\n", iter_wu(i,j),res_wu(i,j))
        end
    end
    time = toc;
    dist_v_b(j,:) = mean(dist(:,:,j),1);
    dist_v_b_wu(j,:) = mean(dist_wu(:,:,j),1);
    dist_v_b_wu2(j,:) = mean(dist_wu2(:,:,j),1);

    fprintf("t=%.2fsec for b=%.2f, Avg TV=%f, Avg KL=%f\n",time, b_space(j), dist_v_b(j,1), dist_v_b(j,2));
end
filename = sprintf('diag_bias_seed%d_%s.mat', seed, datestr(now,'HH_MM_SS_FFF'));
save(filename);
function dist = tv_tringle(s1,b1,s2,b2)
    dist = sum(abs(normcdf(-b1./s1) - normcdf(-b2./s2)));
    diff = @(x) abs(normpdf(x, b1, s1) - normpdf(x, b2, s2));
    dist  = dist + sum(integral(diff, 0, max([b1;b2]) + 6*max([s1;s2]), 'ArrayValued', true));
end
end