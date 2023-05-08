function [sigma_hat, b_hat] = main(samples)
% samples: dxn vector generated from D(W^*, b^*) for some non-negative b^*
[d, n] = size(samples);
sigma_hat = zeros(d,d);
b_hat = zeros(d,1);
% estimate the row norms of W^* and b^*
for i = 1:d
    S = samples(i,:);
    [b_hat(i), sigma_hat(i,i)] = NormBiasEst(S(S>0));
    b_hat(i) = max(0, b_hat(i));
end
% estimate the angles between any two row vectors of W^*
for i = 1:d
    for j = (i+1):d
        theta_hat = pi-2*pi*sum(double(samples(i,:)>b_hat(i)).*double(samples(j,:)>b_hat(j)))/n;
        sigma_hat(i,j) = sqrt(sigma_hat(i,i)*sigma_hat(j,j))*cos(theta_hat);
        sigma_hat(j,i) = sigma_hat(i,j);
    end
end
    