function [mu_hat, sigma2_hat] = NormBiasEst(samples)
% samples: 1xn vector, samples from N(mu^*, sigma2^*, >0)
n = size(samples,2);
% shift and rescale the samples
mu_0 = sum(samples)/n;
sigma2_0 = sum((samples-mu_0).^2)/n;
samples = (samples-mu_0)./sqrt(sigma2_0);
% split the samples into B batches
B = 1;
batch_size = floor(n/B);
v = zeros(2,B);
for b = 1:B
    batch = samples((batch_size*(b-1)+1):batch_size*b);
    v(:,b) = projSGD(batch, -mu_0/sqrt(sigma2_0));
end
best_v = v(:,1);
for b = 1:B
    if b == 1
        best_dist = norm(v-repmat(v(:,b),1,B), 'fro');
    else
        dist = norm(v-repmat(v(:,b),1,B), 'fro');
        if dist < best_dist
            best_dist = dist;
            best_v = v(:,b);
        end
    end
end
sigma2_hat = sigma2_0/best_v(1,1);
mu_hat = mu_0 + best_v(2,1)*sqrt(sigma2_0)/best_v(1,1);
