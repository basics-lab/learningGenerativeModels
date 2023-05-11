%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes an approximation to the gradient of a trucated 
% non-negative normal distribution in the natural parameters U (varaince)
% and V = U*mu, where mu is the mean
%
% grad_lik_function(Y,U,V)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gradU, gradv] = grad_lik_function(Y,U,V)
    [~,n] = size(Y);
    S = Y <= 0;
    Sc =~S;
    try L = chol(U);
    catch
        keyboard;
    end
    invL = inv(L);
    invU = invL*invL';
    gradU = n*(0.5*invU + 0.25*(invU*(V*V')*invU));
    gradv = -n*0.5*invU*V;
    for i=1:n
        si = S(:,i);
        d_si = nnz(si);
        sci = Sc(:,i);
        y_sc = Y(sci, i);
        gradU(sci,sci) = gradU(sci,sci) + -y_sc*y_sc';
        gradv(sci) = gradv(sci) - y_sc;
        if d_si == 0
            continue
        else
            Us = U(si,si);
            invM = inv(chol(Us));
            invUs = invM*invM';
            mu = -invUs*(U(si, sci)*y_sc + 0.5*V(si));
            x = sampleMultiTrunGauss(mu,invUs/2, n);
            % Gradient Updates
            mu_x = mean(x,2);
            gradU(si,si) = gradU(si,si) + -x*x'/n;
            grad_sc_s = -y_sc*mu_x';
            gradU(sci,si) = gradU(sci,si) + grad_sc_s;
            gradU(si,sci) = gradU(si,sci) + grad_sc_s';
            gradv(si) = gradv(si) - mu_x;
        end
    end
end