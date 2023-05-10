clear;clc;
rng(0)
d = 10;
n = 2;
A = eye(d) + randn(d,d)./(2*d);
sigma = A*A';
mu = zeros(d,1);
Y = max(0, mvnrnd(mu,sigma,n)');
U = inv(sigma)/2;
V = -2*U*mu;
grad = grad_lik_function(Y,U,V);
