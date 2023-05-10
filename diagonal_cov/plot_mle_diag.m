%% Setup
clear;clc;
rng(0)
d = 3;
N= 10; % number of points
b_space = linspace(1e-3,4,N); %grid
n_mc = 500;
n=5000;
sigma = ones(3,1);
verbose = true;
%%
dist_v_b = zeros(N,2);
dist_v_b_wu = zeros(N,2);
dist = zeros(n_mc,2,N);
dist_wu = zeros(n_mc,2,N);
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
        [u_hat, v_hat, res, iter] = projGD_MLE_diag(Y, u0, v0);
        dist(i,:,j) = dist_prod_distr(u_hat, v_hat, u_star, v_star);
        if b(1) > 0
            [sigma_hat_wu, b_hat_wu] = main_diag(Y, u0./v0, 1./v0);
            v_hat_wu = diag(sigma_hat_wu).^(0.5);
            u_hat_wu = b_hat_wu./v_hat_wu;
            dist_wu(i,:,j) = dist_prod_distr(u_hat_wu, v_hat_wu, u_star, v_star);
        end
        if verbose
            fprintf("Opt took t=%i iterations with residual=%f\n",iter,res)
        end
    end
    time = toc;
    dist_v_b(j,:) = mean(dist(:,:,j),1);
    dist_v_b_wu(j,:) = mean(dist_wu(:,:,j),1);
    fprintf("t=%.2fsec for b=%.2f, Avg TV=%f, Avg KL=%f\n",time, b_space(j), dist_v_b(j,1), dist_v_b(j,2));
end
figure;
plot(b_space, dist_v_b(:,1));
hold on 
plot(b_space, dist_v_b_wu(:,1));
function dist = dist_prod_distr(u1,v1,u2,v2)
    dist(1) = compute_TV(v1.^(-1), u1./v1, v2.^(-1), u2./v2);
    dist(2) = (0.5*KL_div(u1./v1, diag(v1.^(-2)), u2./v2, diag(v2.^(-2)))).^(0.5);
    %dist(3) = tv_tringle(v1.^(-1), u1./v1, v2.^(-1), u2./v2);
end

function dist = tv_tringle(s1,b1,s2,b2)
    dist = sum(abs(normcdf(-b1./s1) - normcdf(-b2./s2)));
    diff = @(x) abs(normpdf(x, b1, s1) - normpdf(x, b2, s2));
    dist  = dist + sum(integral(diff, 0, max([b1;b2]) + 6*max([s1;s2]), 'ArrayValued', true));
end