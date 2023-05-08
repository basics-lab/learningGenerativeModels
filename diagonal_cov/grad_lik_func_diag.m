function grad = grad_lik_func_diag(Y,u,v)
    [d,n] = size(Y);
    S = Y > 0;
    Sc = Y <= 0;
    %% Gradient WRT u
    gradu = zeros(d,n);
    gradv = zeros(d,n);
    for i =1:d
        y = Y(i, S(i,:));
        gradu(i, S(i,:)) = v(i)*y - u(i);
        gradu(i, Sc(i,:)) = -normpdf(-u(i))/normcdf(-u(i));
        gradv(i, S(i,:)) = 1/v(i) - (v(i)*y - u(i)).*y;
    end
    grad = sum([gradu;gradv],2);
end