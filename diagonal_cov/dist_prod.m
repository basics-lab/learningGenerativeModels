function dist = dist_prod(u1,v1,u2,v2)
    dist(1) = compute_TV(v1.^(-1), u1./v1, v2.^(-1), u2./v2);
    dist(2) = (0.5*KL_div(u1./v1, diag(v1.^(-2)), u2./v2, diag(v2.^(-2)))).^(0.5);
    dist(3) = (0.5*KL_div(u2./v2, diag(v2.^(-2)), u1./v1, diag(v1.^(-2)))).^(0.5);
    dist(4) = (0.5*KL_diag(u2,v2,u1,v1))^0.5;
    dist(5) = (0.5*KL_diag(u1,v1,u2,v2))^0.5;
    dist(6) = tv_tringle(v1.^(-1), u1./v1, v2.^(-1), u2./v2);
    function dist = tv_tringle(s1,b1,s2,b2)
    dist = sum(abs(normcdf(-b1./s1) - normcdf(-b2./s2)));
    diff = @(x) abs(normpdf(x, b1, s1) - normpdf(x, b2, s2));
    dist  = dist + sum(integral(diff, 0, max([b1;b2]) + 6*max([s1;s2]), 'ArrayValued', true));
    end
end