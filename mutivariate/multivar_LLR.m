function lik = multivar_LLR(Y, U1, V1,U2,V2)
    options = statset('TolFun',1e-8);
    [d,n] = size(Y);
    S = Y <= 0;
    Sc =~S;
    L = chol(U1);
    logdetU1 = 2*sum(log(diag(L)));
    L = chol(U2);
    logdetU2 = 2*sum(log(diag(L)));
    lik = zeros(n,1);
    for i=1:n
        s = S(:,i);
        d_s= nnz(S(:,i));
        sc = Sc(:,i);
        lik(i) = -Y(sc,i)'*V1(sc) -Y(sc,i)'*U1(sc,sc)*Y(sc,i) ...
              - (-Y(sc,i)'*V2(sc) -Y(sc,i)'*U2(sc,sc)*Y(sc,i));
        if d_s == 0
            continue
        else
            r1 = U1(s, sc)*Y(sc,i) + 0.5*V1(s);
            r2 = U2(s, sc)*Y(sc,i) + 0.5*V2(s);
            U1s = U1(s,s);
            M1 = chol(U1s);
            invM1 = inv(M1);
            invU1s = invM1*invM1';
            logdetU1s = 2*sum(log(diag(M1)));
            U2s = U2(s,s);
            M2 = chol(U2s);
            invM2 = inv(M2);
            invU2s = invM2*invM2';
            logdetU2s = 2*sum(log(diag(M2)));
            mu1 = -invU1s*r1;
            mu2 = -invU2s*r2;
            lik(i) = lik(i) -r1'*mu1 - 0.5*logdetU1s ...
                          - (-r2'*mu2 - 0.5*logdetU2s) ...
            + log(mvncdf(zeros(d_s,1),mu1,0.5*invU1s,options)/mvncdf(zeros(d_s,1),mu2,0.5*invU2s,options));
        end
    end
    lik = n*(-0.25*V1'*(U1\V1) + 0.5*logdetU1) - n*(-0.25*V2'*(U2\V2) + 0.5*logdetU2) ...
          + sum(lik);
end