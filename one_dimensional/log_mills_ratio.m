function r = log_mills_ratio(x)
    S = x < 20;
    Sc = ~S;
    r(S)  = log(normcdf(-x(S))) + 0.5*sqrt(2*pi) + 0.5*x(S).^2;
    r1 = 1./x(Sc);
    r(Sc) = r1;
    n=30;
    for k=0:n
        r1 = (-1)*(2*k + 1)*r1.*x(Sc).^(-2);
        r(Sc) = r(Sc) + r1;
    end
    r(Sc) = log(r(Sc));
end