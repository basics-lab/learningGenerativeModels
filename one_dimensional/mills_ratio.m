function r = mills_ratio(x)
    S = x < 20;
    Sc = ~S;
    r(S)  = normcdf(-x(S))./normpdf(x(S));
    r1 = 1./x(Sc);
    r(Sc) = r1;
    n=20;
    for k=0:n
        r1 = (-1)*(2*k + 1)*r1.*x(Sc).^(-2);
        r(Sc) = r(Sc) + r1;
    end
end