function run_diag_dim(N, n, n_mc, seed)
    rng(seed);
    if getenv('USER') == "justinkang"
        parpool(str2num(getenv('SLURM_CPUS_ON_NODE')));
        fprintf("ON SLURM, CREATING A BIGGER POOL\n")
    end
    KL_arr = zeros(N,n_mc);
    d_space = ceil(logspace(log10(5),log10(20),N)); %number of samples
    filename = sprintf('diag_dim_seed%d_%s.mat', seed, datestr(now,'HH_MM_SS_FFF'));
    for j=1:N
        d = d_space(j);
        fprintf("Running now with d=%i\n", d);
        sigma = ones(d,1);
        b = ones(d,1);
        u_star = b./sigma;
        v_star = 1./sigma;
        for i=1:n_mc
            % Param
            u0 = b + randn(d,1);
            v0 = 0.5 + exprnd(0.5, d,1);
            % Sampling
            Y = max(0, ((sigma*ones(1,n)).*randn(d,n)) + b);
            % Computing Estimates
            [u_hat,v_hat,res,iter] = projGD_MLE_diag(Y, u0, v0);
            fprintf("Computed Optimal in t=%.3f iter with res =%f\n", iter, res)
            %[sigma_hat_wu, b_hat_wu,res, iter, ~] = main_diag(Y, u0./v0,1./v0);
            %fprintf("Computed Optimal Wu in t=%.3f iter with res =%f\n", iter, res)
            %v_hat_wu = diag(sigma_hat_wu).^(-0.5);
            %u_hat_wu = b_hat_wu'.*v_hat_wu;
            KL_arr(j,i) = KL_diag(u_hat,v_hat, u_star, v_star);
        end
        save(filename);
    end
end