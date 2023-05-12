function [u, v] = cvx_solve(X, y)
    [d,~] = size(X);
    S = y > 0;
    Sc = y <= 0;
    absS = nnz(S);
    cvx_begin
        variable v
        variable u(d)
         maximize(-0.5*sum_square(v*y(S) - X(:,S)'*u) + ...
             absS*log(v) + sum(log_normcdf(-u'*X(:,Sc))));
    cvx_end
end
