%% Univariate likelihood function
function l = lik_func( u, v, X, y)
 S = y > 0;
 Sc = y <= 0;
 l = sum(-0.5*(v*y(S)' - u'*X(:,S)).^2 + log(v)) + ...
     sum(log(normcdf(-u'*X(:,Sc))));
end

function l = lik_func( u, v, X, y) 
% Precompute some values to avoid repeated calculations 
z = u'X; 
vy = vy; 
% Use element-wise operations and vectorization to speed up computation 
l = sum(-0.5*(vy(z > 0) - z(z > 0)).^2 + log(v)) + ...
    sum(log(normcdf(-z(z <= 0))));
end