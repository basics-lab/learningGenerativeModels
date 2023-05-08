function l = lik_func_diag(Y, u, v)
    [d, n] = size(Y);
    S = Y > 0;
    rowsumS = sum(S,2);
    rowsumSc = n - rowsumS;
    dist = Y .* v(:, ones(n,1)) - u;
    l = rowsumS'*log(v) + rowsumSc'*log(normcdf(-u)) + -0.5*sum(dist(S).^2);
    if isnan(l)
        l = -Inf;
    end
end
