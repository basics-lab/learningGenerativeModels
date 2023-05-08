function v_bar = projSGD(samples, thres)
% samples: 1xn vector, samples from N(mu^*, sigma2^*, >thres)
lambda = 0.1;
v = [1;0];
T = size(samples,2);
r = 3;
v_avg = [0;0];
batch_size = 10;
num_batch = floor(T/batch_size);
for t = 1:num_batch
    mu = v(2,1)/v(1,1);
    sigma2 = 1/v(1,1);
    z = sampleTrunGauss(mu, sigma2, thres, batch_size);
    x = samples((batch_size*(t-1)+1):batch_size*t);
    grad = [0.5*mean(x.^2) - 0.5*mean(z.^2); mean(z)-mean(x)];
    v = v - grad/(lambda*t);
    if v(1,1) > r
        v(1,1) = r;
    elseif v(1,1) < 1/r
        v(1,1) = 1/r;
    end
    if v(2,1) > r
        v(2,1) = r;
    elseif v(2,1) < -r
        v(2,1) = -r;
    end
    v_avg = v_avg + v;
end
v_bar = v_avg/num_batch;