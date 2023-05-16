function run_condition_number(N, n, n_mc, seed)
    rng(seed);
    if getenv('USER') == "justinkang"
        parpool(str2num(getenv('SLURM_CPUS_ON_NODE')));
        fprintf("ON SLURM, CREATING A BIGGER POOL\n")
    end
    d = 3;
    kappa = linspace(1, 25^2,N);
    b = ones(d,1);
    dist_v_cond = zeros(n_mc, 6, N);
    dist_v_cond_wu = zeros(n_mc, 6, N);
    filename = sprintf('diag_cond_seed%d_%s.mat', seed, datestr(now,'HH_MM_SS_FFF'));
    for j=1:N
        fprintf("Running now with kappa=%.2f\n", kappa(j));
        sigma2_max = sqrt(kappa(j));
        sigma2_min = sqrt(kappa(j))^(-1);
        keyboard;
        for i=1:n_mc
            % Param
            sigma = [sigma2_max;sigma2_min;1].^0.5;
            u_star = b./sigma;
            v_star = 1./sigma;
            u0 = b + randn(d,1);
            v0 = 0.5 + exprnd(0.5, d,1);
            % Sampling
            Y = max(0, ((sigma*ones(1,n)).*randn(d,n)) + b);
            % Computing Estimates
            [u_hat,v_hat,res,iter] = projGD_MLE_diag(Y, u0, v0);
            fprintf("Computed Optimal in t=%.3f iter with res =%f\n", iter, res)
            [sigma_hat_wu, b_hat_wu,res, iter, ~] = main_diag(Y, u0./v0,1./v0);
            fprintf("Computed Optimal Wu in t=%.3f iter with res =%f\n", iter, res)
            v_hat_wu = diag(sigma_hat_wu).^(-0.5);
            u_hat_wu = b_hat_wu'.*v_hat_wu;
            dist_v_cond(i,:,j) = dist_prod(u_hat,v_hat, u_star,v_star);
            dist_v_cond_wu(i,:,j) = dist_prod(u_hat_wu, v_hat_wu, u_star,v_star);
        end
        save(filename);
    end
end
