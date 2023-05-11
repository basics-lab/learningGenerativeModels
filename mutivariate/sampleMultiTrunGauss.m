function samples = sampleMultiTrunGauss(mu, sigma2, num_samples)
% sample from a truncated normal distribution N(mu, sigma2, >thres)
d = length(mu);
samples = [];
T = max(100, num_samples*10);
while size(samples,2)<num_samples
    x = mvnrnd(mu,sigma2, T)';
    samples = [samples, x(:,sum(x < 0,1) == d)];
end
if num_samples == 1
    samples = samples(:,1);
else
    samples = samples(:,1:num_samples);
end