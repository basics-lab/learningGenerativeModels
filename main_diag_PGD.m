function [sigma_hat, b_hat, zero_sample] = main_diag_PGD(samples)
% samples: dxn vector generated from D(W^*, b^*) for some non-negative b^*
[d, n] = size(samples);
zero_sample = false;
%sigma_hat = diag(s0.^2);
%b_hat = b0;
% estimate the row norms of W^* and b^*
for i = 1:d
    S = samples(i,:);
    untrunc = S(S>0);
    zero_sample = zero_sample | isempty(untrunc);
    [b_hat(i), sigma_hat(i,i)] = NormBiasEst(untrunc);
end