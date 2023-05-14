function [samples, failed] = sampleMultiTrunGauss(mu, sigma2, num_samples)
% sample from a truncated normal distribution N(mu, sigma2, >thres)
d = length(mu);
samples = [];
T = max(100, num_samples*10);
t = 0;
failed=false;
while size(samples,2)<num_samples
    x = mvnrnd(mu,sigma2, T)';
    samples = [samples, x(:,sum(x < 0,1) == d)];
    t =t+1;
    if(t == 100000)
        fprintf("Gradient took too long, moving on...\n");
        failed = true;
        return
    end
end
if num_samples == 1
    samples = samples(:,1);
else
    samples = samples(:,1:num_samples);
end