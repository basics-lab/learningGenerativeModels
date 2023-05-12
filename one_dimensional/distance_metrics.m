function [dist] = distance_metrics(u1, v1, u2, v2, X)
    [~, n] = size(X);
    w1 = u1/v1; s1 = 1/v1;
    w2 = u2/v2; s2 = 1/v2;
    std = max([s1,s2]);
    mu1 = X'*w1;
    mu2 = X'*w2;
%     % KL Div
%     div = 0.5*(mu1 - mu2).^2/s2^2 + 0.5*s1^2/s2^2 + log(s2/s1) - 0.5;
%     dist(1) = sqrt(0.5*mean(div));
%     % Standard TV
%     diff = @(x) abs(normpdf(x, zeros(n,1), s1) - normpdf(x, mu1-mu2, s2));
%     tv_dist = integral(diff, -6*std, max(mu1-mu2) + 6*std, 'ArrayValued', true);
%     dist(2) =  0.5*mean(tv_dist);
    % Truncated TV
    upper = max([mu1;mu2]) + 6*std;
    diff = @(x) abs(normpdf(x, mu1, s1) - normpdf(x, mu2, s2));
    tv_dist_trunc = integral(diff, 0, upper, 'ArrayValued', true) + ...
           abs(normcdf(-mu1/s1) - normcdf(-mu2/s2));
    dist =  0.5*mean(tv_dist_trunc);
end