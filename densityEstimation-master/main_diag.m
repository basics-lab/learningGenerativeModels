function [sigma_hat, b_hat] = main_diag(samples,b0,s0)
% samples: dxn vector generated from D(W^*, b^*) for some non-negative b^*
[d, n] = size(samples);
sigma_hat = diag(s0.^2);
b_hat = b0;
% estimate the row norms of W^* and b^*
for i = 1:d
    S = samples(i,:);
    [b_hat(i), sigma_hat(i,i)] = NormBiasEst(S(S>0));
    b_hat(i) = max(0, b_hat(i));
end