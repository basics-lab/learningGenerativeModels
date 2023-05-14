function KL = KL_diag(u1,v1,u2,v2)
    d = length(u1);
    n=1000000;
    sigma = v2.^(-1);
    b = u2./v2;
    Y = max(0, ((sigma*ones(1,n)).*randn(d,n)) + b);
    KL = -n^(-1)*LLR_diag(Y,u1,v1,u2,v2);
    if KL < 0
        fprintf("[Error] KL is negative, probably we need to take more samples\n")
    end
end