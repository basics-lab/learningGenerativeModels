function dist = dist_prod_distr(u1,v1,u2,v2)
    dist(1) = compute_TV(v1.^(-1), u1./v1, v2.^(-1), u2./v2);
    dist(2) = (0.5*KL_div(u1./v1, diag(v1.^(-2)), u2./v2, diag(v2.^(-2)))).^(0.5);
    %dist(3) = tv_tringle(v1.^(-1), u1./v1, v2.^(-1), u2./v2);
end