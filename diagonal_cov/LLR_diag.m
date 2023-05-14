function l = LLR_diag(Y, u1, v1,u2,v2)
    [d, n] = size(Y);
    S = Y > 0;
    rowsumS = sum(S,2);
    rowsumSc = n - rowsumS;
    dist1 = Y .* v1(:, ones(n,1)) - u1;
    dist2 = Y .* v2(:, ones(n,1)) - u2;
    l = rowsumS'*(log(v1./v2)) + rowsumSc'*log(normcdf(-u1)./normcdf(-u2)) + -0.5*sum(dist1(S).^2) + 0.5*sum(dist2(S).^2);
    if isnan(l)
        l = -Inf;
    end
end
