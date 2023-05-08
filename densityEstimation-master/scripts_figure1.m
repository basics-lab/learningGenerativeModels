%-----------subfigure 1: err vs #samples---------------%
num_samples = 20000:30000:200000;
d = 5;
k = 5;
num_runs = 10;
ns = size(num_samples, 2);
err_b = zeros(ns, num_runs);
err_sigma = zeros(ns, num_runs);
sqrt_kl = zeros(ns, num_runs); % upper bound on TV distance
for i = 1:ns
    num_sample = num_samples(1,i)
    for j = 1:num_runs
        b_star = max(0, randn(d,1));
        W_star = randn(d,k)/sqrt(k);
        [U,~,~] = svd(W_star);
        W_star = U;
        Z = randn(k, num_sample);
        samples = max(0, W_star*Z + b_star);
        [sigma_hat, b_hat] = main(samples);
        sigma_star = W_star*W_star';
        err_sigma(i, j) = norm(sigma_hat - sigma_star, 'fro')/norm(W_star,'fro')^2;
        err_b(i, j) = norm(b_hat - b_star)/norm(W_star, 'fro');
        kl_div_1 = KL_div(b_star, sigma_star, b_hat, sigma_hat); % can be complex if sigma_hat/sigma_star is not PSD
        kl_div_2 = KL_div(b_hat, sigma_hat, b_star, sigma_star);
        sqrt_kl(i, j) = sqrt(min(kl_div_1, kl_div_2)/2); % upper bound on TV distance
    end
end
save('fig1_subfig1.mat','num_samples', 'err_b', 'err_sigma', 'sqrt_kl')
%-----------subfigure 2: err vs dimension---------------%
clear all;
num_sample = 500000;
dims = 5:4:25;
num_runs = 10;
num_dims = size(dims,2);
err_b = zeros(num_dims, num_runs);
err_sigma = zeros(num_dims, num_runs);
sqrt_kl = zeros(num_dims, num_runs);
for i = 1:num_dims
    d = dims(1,i)
    for j = 1:num_runs
        b_star = max(0, randn(d,1));
        W_star = randn(d, d)/sqrt(d);
        [U,~,~] = svd(W_star);
        W_star = U; % condition number = 1 
        Z = randn(d, num_sample);
        samples = max(0, W_star*Z + b_star);
        [sigma_hat, b_hat] = main(samples);
        sigma_star = W_star*W_star';
        err_sigma(i, j) = norm(sigma_hat - sigma_star, 'fro')/norm(W_star,'fro')^2;
        err_b(i, j) = norm(b_hat - b_star)/norm(W_star, 'fro');
        kl_div_1 = KL_div(b_star, sigma_star, b_hat, sigma_hat); % can be complex if sigma_hat/sigma_star is not PSD
        kl_div_2 = KL_div(b_hat, sigma_hat, b_star, sigma_star);
        sqrt_kl(i, j) = sqrt(min(kl_div_1, kl_div_2)/2); % upper bound on TV distance
    end
end
save('fig1_subfig2.mat', 'dims', 'err_b', 'err_sigma', 'sqrt_kl')
%-----------subfigure 3: err vs condition number---------------%
clear all;
num_sample = 500000;
kappas = 1:4:25; % condition number of WW'
d = 5;
num_runs = 10;
num_kappas = size(kappas,2);
err_b = zeros(num_kappas, num_runs);
err_sigma = zeros(num_kappas, num_runs);
sqrt_kl = zeros(num_kappas, num_runs);
for i = 1:num_kappas
    kappa = kappas(1, i)
    for j = 1:num_runs
        b_star = max(0, randn(d,1));
        W_star = randn(d, d)/sqrt(d);
        [U,~,~] = svd(W_star);
        new_S = [];
        lambda = (sqrt(kappa))^(1/(d-1));
        for k = 1:d
             new_S = [new_S, 1/lambda^(k-1)];
        end
        W_star = U*diag(new_S);
        Z = randn(d, num_sample);
        samples = max(0, W_star*Z + b_star);
        [sigma_hat, b_hat] = main(samples);
        sigma_star = W_star*W_star';
        err_sigma(i, j) = norm(sigma_hat - sigma_star, 'fro')/norm(W_star,'fro')^2;
        err_b(i, j) = norm(b_hat - b_star)/norm(W_star, 'fro');
        kl_div_1 = KL_div(b_star, sigma_star, b_hat, sigma_hat); % can be complex if sigma_hat/sigma_star is not PSD
        kl_div_2 = KL_div(b_hat, sigma_hat, b_star, sigma_star);
        sqrt_kl(i, j) = sqrt(min(kl_div_1, kl_div_2)/2); % upper bound on TV distance
    end
end
save('fig1_subfig3.mat', 'kappas', 'err_b', 'err_sigma', 'sqrt_kl')