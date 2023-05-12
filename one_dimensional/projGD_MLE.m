function [u,v, res, iter] = projGD_MLE(y, X, u0, v0)
    %% Parameters
    max_iter = 10000;
    tol = 0.01;
    tolv = 1e-2;
    beta = 0.4;
    %% Definitions
    S = y > 0;
    Sc = y <= 0;
    gradu = @(u,v) X(:,S) * (v*y(S)' - u'*X(:,S))' ...
                 - X(:,Sc)* (normpdf(-u'*X(:,Sc))./normcdf(-u'*X(:,Sc)))';
    gradv = @(u,v) nnz(S)/v - y(S)'*(v*y(S) - X(:,S)'*u);
    lik = @(u,v,du,dv,a) -lik_func((u + a*du), max([v + a*dv, tolv]), X, y);
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
        u = u + alpha*du;
        v = max([v + alpha*dv, tolv]);
    end
    iter = t;
end