function [mu_hat, sigma2_hat, res, iter] = gradNormWu(y,u0,v0)
    %% Parameters
    max_iter = 10000;
    n = length(y);
    tol = 0.001;
    tolv = 1e-2;
    u_bnd = 8;
    beta = 0.5;
    %% Definitions
    gradu = @(u,v) sum((v*y - u) - normpdf(-u)./(1 - normcdf(-u)));
    gradv = @(u,v) n/v - sum(y.*(v*y - u));
    lik = @(u,v,du,dv,a) -one_dim_trunc_lik(y,min([max([u + a*du, -u_bnd]), u_bnd]),max([v + a*dv, tolv]));
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
        res = (du'*du + dv'*dv);
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
        u = min([max([u + alpha*du, -u_bnd]), u_bnd]);
        v = max([v + alpha*dv, tolv]);
    end
    iter = t;
    function out = one_dim_trunc_lik(y,u,v)
        out = sum(log(v) -0.5*(v*y - u).^2 - log(1 - normcdf(-u)));
         if isnan(out)
            out = -Inf;
        end
    end
sigma2_hat = 1/v.^2;
mu_hat = u/v;
end