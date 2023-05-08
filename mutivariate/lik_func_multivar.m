function lik_func_multivar(mu, V, Y)
    [d, n] = size(Y);
    S = Y <= 0;
    P = V'*V;
    for i =1:n
        s = S(:,i);
        sc = S(:,i) == 0;
        Ps = P(s,s);
        Psc = P(sc,sc);
        Pssc = P(s,sc);
        Pscs = P(sc,s);
    end
    
end