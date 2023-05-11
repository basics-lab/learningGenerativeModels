function lik = lik_func_multivar(Y, U, V)
    options = statset('TolFun',1e-5);
    [d,n] = size(Y);
    S = Y <= 0;
    Sc =~S;
    L = chol(U);
    logdetU = 2*sum(log(diag(L)));
    lik = zeros(n,1);
    for i=1:n
        s = S(:,i);
        d_s= nnz(S(:,i));
        sc = Sc(:,i);
        lik(i) = -Y(sc,i)'*V(sc) -Y(sc,i)'*U(sc,sc)*Y(sc,i) + (d - d_s)*log(2)/2;
        if d_s == 0
            continue
        else
            r = U(s, sc)*Y(sc,i) + 0.5*V(s);
            Us = U(s,s);
            M = chol(Us);
            invM = inv(M);
            invUs = invM*invM';
            logdetUs = 2*sum(log(diag(M)));
            mu = -invUs*r;
            lik(i) = lik(i) -r'*mu - 0.5*logdetUs ...
            + log(mvncdf(zeros(d_s,1),mu,0.5*invUs,options));
        end
    end
    lik = n*(-0.25*V'*(U\V) + 0.5*logdetU) + sum(lik);
end