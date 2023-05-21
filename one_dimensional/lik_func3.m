%% Univariate likelihood function
function l = lik_func3( u, v, X, y)
 S = y > 0;
 Sc = y <= 0;
 l = sum(-0.5*(v*y(S)' - u'*X(:,S)).^2 + log(v)) + ...
     sum(lognormcdf_mills(u'*X(:,Sc)));
 if isnan(l)
    l = -Inf;
 end
    function ret  = lognormcdf_mills(x)
        ret  = log_mills_ratio(x) - 0.5*log(2*pi) - 0.5*x.^2;
    end
end
