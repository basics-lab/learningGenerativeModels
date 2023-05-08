function [u,v,res,iter] = projGD_MLE_diag(Y, u0, v0)
    %% Parameters
    [d,~] = size(Y);
    max_iter = 5000;
    tol = 0.001;
    tolv = 1e-3;
    beta = 0.5;
    grad = @(u,v) grad_lik_func_diag(Y,u,v);
    lik = @(u,v,du,dv,a) -lik_func_diag(Y, (u + a*du), max([v + a*dv, tolv*ones(d,1)],[],2));
    u = u0;
    v = v0;
    t = 0;
    res = inf;
    curr = lik(u,v,0,0,0);
    %% Iterations
    while t < max_iter && (res > tol)
        t = t + 1;
        alpha = 1;
        g = grad(u,v);
        du = g(1:d);
        dv = g(d+1:end);
        % Line Search
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
         u = u + alpha*du;
         v = max([v + alpha*dv, tolv*ones(d,1)], [], 2);
    end
    iter = t;
end