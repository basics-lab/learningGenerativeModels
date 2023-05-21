%% Load Data
clear;clc;
X ={};
Y = {};
myFiles_n1 = dir(fullfile("data_may_16","*_n1000_*.mat"));
myFiles_n2 = dir(fullfile("data_may_16","*_n10000_*.mat"));
n_runs = length(myFiles_n1);
for k = 1:length(myFiles_n1)
    filename_n1 = strcat(myFiles_n1(k).folder,'/',myFiles_n1(k).name);
    filename_n2 = strcat(myFiles_n2(k).folder,'/',myFiles_n2(k).name);
    X{k} = load(filename_n1);
    Y{k} = load(filename_n2);
end

%% Compute Means
[n_pts,n_mc]  = size(X{1}.dist_v_n);
d_space = X{1}.d_space;
dist_v_n1 = zeros(n_pts,1);
dist_v_n2 = zeros(n_pts,1);
for k=1:n_runs
    dist_v_n1 = dist_v_n1 + mean(X{k}.dist_v_n, 2, "omitnan");
    dist_v_n2 = dist_v_n2 + mean(Y{k}.dist_v_n, 2, "omitnan");
end
n1_mean = dist_v_n1/n_runs;
n2_mean = dist_v_n2/n_runs;
%% Compute std
dist_v_n1 = zeros(n_pts,1);
dist_v_n2 = zeros(n_pts,1);
for k=1:n_runs
    dist_v_n1 = dist_v_n1 + mean((X{k}.dist_v_n - n1_mean).^2, 2, "omitnan");
    dist_v_n2 = dist_v_n2 + mean((Y{k}.dist_v_n - n2_mean).^2, 2, "omitnan");
end
n1_var = dist_v_n1/n_runs;
n2_var = dist_v_n2/n_runs;
n1std = n1_var.^0.5;
n2std = n2_var.^0.5;
n1stderr = n1std/(n_mc*n_runs)^0.5;
n2stderr = n2std/(n_mc*n_runs)^0.5;
figure;
a = axes;
errorbar(d_space, n1_mean, n1stderr, "Linewidth", 2);
hold on
errorbar(d_space, n2_mean, n2stderr, "-.", "Linewidth", 2);
a.XScale = "log";
a.YScale = "log";
xlim([6, 80]);
ylim("auto");
xlabel("Input dimension $k$", "Interpreter","latex", "FontSize",16)
ylabel("$d\left((\hat{w},\hat{\sigma}^2), ( w^*,\sigma^{2*}) \right)$", "Interpreter","latex", "FontSize",16)
legend("$n = 10^3$","$n = 10^4$", "Location", "best", "Interpreter","Latex","FontSize",16)