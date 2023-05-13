function dist = dist_prod_distr(u1,v1,u2,v2)
    dist(1) = compute_TV(v1.^(-1), u1./v1, v2.^(-1), u2./v2);
    dist(2) = (0.5*KL_div(u1./v1, diag(v1.^(-2)), u2./v2, diag(v2.^(-2)))).^(0.5);
    dist(3) = norm(u1./v1 - u2./v2);
    dist(4) = norm(v1.^(-2) - v2.^(-2));
    dist(5) = norm(u1 - u2);
    dist(6) = norm(v1 - v2);
    dist(7) = norm(u1.*v1 - u2.*v2);
    dist(8) = norm(v1.^2 - v2.^2);
    %dist(3) = tv_tringle(v1.^(-1), u1./v1, v2.^(-1), u2./v2);
end