function [u,v] = projGD_MLE(y, X, u0, v0)
    %% Parameters
    max_iter = 100000;
    step_max = 1e-3;
    step_min = 9e-6;
    tol = 1e-2;
    %% Definitions
    S = y > 0;
    Sc = y <= 0;
    options = optimset('TolX', step_min);
    gradu = @(u,v) X(:,S) * (v*y(S)' - u'*X(:,S))' ...
                 - X(:,Sc)* (normpdf(-u'*X(:,Sc))./normcdf(-u'*X(:,Sc)))';
    gradv = @(u,v) nnz(S)/v - y(S)'*(v*y(S) - X(:,S)'*u);
    lik = @(u,v,du,dv,a) -lik_func((u + a*du), (v + a*dv), X, y);
    u = u0;
    v = v0;
    t = 0;
    du = Inf;
    dv = Inf;
    %% Iterations
    while t < max_iter && (norm(du) + abs(dv) > tol)
        t = t + 1;
        du = gradu(u,v);
        dv = gradv(u,v);
        step = fminsearch(@(a) lik(u,v,du,dv,a) + 1e10*(a < 0), 0, options); % Line Search
        step = max([step, step_min]);
        u = u + (step)*du;
        v = max([v + (step)*dv, 1e-3]);
    end
end