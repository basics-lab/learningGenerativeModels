function [mu_hat, sigma2_hat, res, iter] = gradNormWu2(y,u0,v0)
    %% Parameters
    max_iter = 10000;
    n = length(y);
    tol = 0.01;
    tolv = 1e-2;
    u_upper = 8;
    u_lower = -8;
    beta = 0.5;
    %% Definitions
    gradu = @(u,v) sum((v*y - u) - v^(-0.5)*normpdf(-u/v^0.5)./(1 - normcdf(-u/v^0.5)));
    gradv = @(u,v) n/v - sum(y.*(v*y - u)) + 0.5*u/v^(-1.5)*normpdf(-u/v^0.5)./(1 - normcdf(-u/v^0.5));
    lik = @(u,v,du,dv,a) -one_dim_trunc_lik(y,min([max([u + a*du, -8]), 8]),max([v + a*dv, tolv]));
    u = u0;
    v = v0;
    t = 0;
    res = inf;
    curr = lik(u,v,0,0,0);
    %% Iterations
    while t < max_iter && (res > tol)
        t = t + 1;
        alpha = 1;
        du = gradu(u,v);
        dv = gradv(u,v);
        res = du'*du + dv'*dv;
        diff = lik(u,v,du,dv,alpha)-curr;
        if isnan(diff)
            diff = Inf;
        end
        while diff > -alpha/2*res
            alpha = beta*alpha;
            diff = lik(u,v,du,dv,alpha) - curr;
            if isnan(diff)
                diff = Inf;
            end
        end
        curr = diff + curr;
        u = min([max([u + alpha*du, -8]), 8]);
        v = max([v + alpha*dv, tolv]);
    end
    iter = t;
    function out = one_dim_trunc_lik(y,u,v)
        u =
        out = sum(log(v) -0.5*(v*y - u).^2 - log(1 - normcdf(-u)));
         if isnan(out)
            out = -Inf;
        end
    end
sigma2_hat = 1/v.^2;
mu_hat = u/v;
end
