%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes an approximaiton to the gradient of a trucated 
% non-negative normal distribution with 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gradU, gradv] = grad_lik_function(Y,U,V)
    [~,n] = size(Y);
    S = Y <= 0;
    Sc =~S;
    invU = inv(U);
    gradU = n*0.5*invU;
    gradU = gradU + n*0.25*(invU*V*V'*invU);
    gradv = -n*0.5*invU*V;
    n_mc = 100000;
    for i=1:n
        % For readability
        si = S(:,i);
        sci = Sc(:,i);
        y_sc = Y(sci, i);
        Us = U(si,si);
        invM = inv(chol(Us));
        invUs = invM*invM';
        mu = invUs*(U(si, sci)*y_sc + 0.5*V(si));
        x = mvnrnd(mu,invUs/2, n_mc)';
        % Check which samples to reject
        x = x(:,sum(x < 0,1) == nnz(si));
        [~,n_acc] = size(x);
        % Gradient Updates
        gradU(sci,sci) = gradU(sci,sci) + -y_sc*y_sc';
        gradU(si,si) = gradU(si,si) + -x*x';
        grad_sc_s = -repmat(y_sc,1,n_acc)*x';
        gradU(sci,si) = gradU(sci,si) + grad_sc_s;
        gradU(si,sci) = gradU(si,sci) + grad_sc_s';
        gradv(sci) = gradv(sci) - y_sc;
        gradv(si) = -mean(x,2);
    end
end