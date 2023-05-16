%% Univariate likelihood function
function l = lik_func2( u, v, X, y)
 S = y > 0;
 Sc = y <= 0;
 l = sum(-0.5*(v*y(S)' - u'*X(:,S)).^2 + log(v)) + ...
     sum(lognormcdf(-u'*X(:,Sc)));
 if isnan(l)
    l = -Inf;
 end
end