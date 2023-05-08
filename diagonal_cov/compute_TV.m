%%
function l = compute_TV(sigma1,b1,sigma2,b2)
    dom = max([b1;b2]) + 5*max([sigma1;sigma2])*ones(1,3);%bounds(sigma1,b1,sigma2,b2,100);% Integration bounds
    d = 3;
    surfaces = (dec2bin(0:2^d-1)' - '0') == 1;
    l = 0;
    set = 1:3;
    for i=1:length(surfaces)
        idx = surfaces(:,i)';
        if nnz(idx) == 0
            l = l + integrand(0,0,0,b1,sigma1,b2,sigma2, idx);
        elseif nnz(idx) == 1
            if dom(set(idx)) < 1e-3
                continue
            end
            f = @(x) integrand(x,0,0,b1,sigma1,b2,sigma2, idx);
            l = l + integral(f,0,dom(set(idx)));
        elseif nnz(idx) == 2
            s = set(idx);
            if dom(s(1)) < 1e-3 || dom(s(2)) < 1e-3
                continue
            end
            f = @(x,y) integrand(x,y,0,b1,sigma1,b2,sigma2, idx);
            l = l + integral2(f,0,dom(s(1)),0,dom(s(2)));
        elseif nnz(idx) == 3
            if dom(1) < 1e-3 || dom(2) < 1e-3 || dom(3) < 1e-3 
                continue
            end
            f = @(x,y,z) integrand(x,y,z,b1,sigma1,b2,sigma2, [true,true,true]);
            l = l + integral3(f,0,dom(1),0,dom(2),0,dom(3));
        end
    end
    l = 0.5*l;
    %% Bounds 
    function dom = bounds(sigma1, mu1, sigma2, mu2, tol)
        b1 = mu1 + tol*sigma1;
        b2 = mu2 + tol*sigma2;
        dom = max([b1';b2';zeros(1,3)]);
    end
    %% TV Integrand
    function l = integrand(x,y,z,b1,sigma1,b2,sigma2, idx)
        [nrow, ncol] = size(x);
        x = x(:);
        y = y(:);
        z = z(:);
        if nnz(idx) == 3
            l = abs(1/prod(sigma1)*mvnpdf(([x,y,z] - b1')./sigma1') - 1/prod(sigma2)*mvnpdf(([x,y,z]-b2')./sigma2'));
        elseif nnz(idx) == 2
            l = abs(normcdf(-b1(~idx)/sigma1(~idx))/prod(sigma1(idx))*mvnpdf(([x,y] - b1(idx)')./sigma1(idx)') - ...
                    normcdf(-b2(~idx)/sigma2(~idx))/prod(sigma2(idx))*mvnpdf(([x,y] - b2(idx)')./sigma2(idx)'));
            l = reshape(l,[nrow,ncol]);
        elseif nnz(idx) == 1
            l = abs(prod(normcdf(-b1(~idx)./sigma1(~idx)))/sigma1(idx)*normpdf((x - b1(idx)')./sigma1(idx)') - ...
                    prod(normcdf(-b2(~idx)./sigma2(~idx)))/sigma2(idx)*normpdf((x - b2(idx)')./sigma2(idx)'));
        else
            l = abs(prod(normcdf(-b1./sigma1)) - prod(normcdf(-b2./sigma2)));
        end
        l = reshape(l,[nrow,ncol]);
    end
end
%% Potentially useful function for Monte-Carlo (does not work)
% function f = tv_dist(X, sigma1, mu1, sigma2, mu2)
%     S = X <= 0;
%     Sc= ~S;
%     X(S) = 0;
%     X1 = (X - mu1)./sigma1;
%     X2 = (X - mu2)./sigma2;
%     X1(S) = normcdf(X1(S));
%     X2(S) = normcdf(X2(S));
%     X1(Sc) = normcdf(X1(Sc));
%     X2(Sc) = normcdf(X2(Sc));
%     X1 = prod(X1, 2);
%     X2 = prod(X2, 2);
%     f = abs(X1 - X2);
% end