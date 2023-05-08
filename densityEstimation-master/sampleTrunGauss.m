function samples = sampleTrunGauss(mu, sigma2, thres, num_samples)
% sample from a truncated normal distribution N(mu, sigma2, >thres)
samples = [];
T = max(100, num_samples*10);
while size(samples,2)<num_samples
    Z = normrnd(mu, sqrt(sigma2), 1, T);
    samples = [samples, Z(Z>thres)];
end
if num_samples == 1
    samples = samples(1,1);
else
    samples = samples(1,1:num_samples);
end