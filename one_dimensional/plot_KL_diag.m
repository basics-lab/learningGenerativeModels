%% Load files
clear;clc;
myFiles = dir(fullfile("data_diag_cov_dim",'*.mat'));
X ={};
n_runs = length(myFiles);
for k = 1:length(myFiles)
    filename = strcat(myFiles(k).folder,'/',myFiles(k).name);
    X{k} = load(filename);
end
N = X{1}.N;
n_mc = X{1}.n_mc;
d_space = X{1}.d_space;
dist_v_dim = zeros(N,1);
%% Mean
for i=1:n_runs
    dist_v_dim = dist_v_dim + mean(X{i}.KL_arr, 2, "omitnan");
end
mean_dist = dist_v_dim/n_runs;
%% Variance
var_v_dist = zeros(N,1);
for i=1:n_runs
    for j=1:N
        var_v_dist = var_v_dist + mean((X{i}.KL_arr - mean_dist).^2, 2,"omitnan");
    end
end
emp_std = (var_v_dist/n_runs).^0.5;
std_err = emp_std/sqrt(n_mc*n_runs);
figure;
a = axes;
hold on;
errorbar(d_space, mean_dist, std_err, "LineWidth", 2);
a.XScale = "log";
a.YScale = "log";
% x = [kappa_space, fliplr(kappa_space)];
% curve1 = [mean_dist_wu + 0.1*var_v_cond_wu, fliplr(mean_dist_wu - 0.1*var_v_cond_wu)];
% fill(x,curve1(1,:),[0.55, 0.81, 0.89]);
% hold on;
% curve2 = [mean_dist + 0.1*var_v_cond, fliplr(mean_dist - 0.1*var_v_cond)];
% fill(x,curve2(1,:),[1, 0.4, 0.2]);
%plot(kappa_space, mean_dist(1,:));
%plot(kappa_space, mean_dist_wu(1,:));
%xlim([1 24^2])
%ylim([0.012 0.032])
xlabel("Condition Number of $\Sigma^*$", "Interpreter","latex", "FontSize",16)
ylabel("$d\left((\hat{w},\hat{\Sigma}), ( w^*,\Sigma^*) \right)$", "Interpreter","latex", "FontSize",16)
legend("Our TV", "Location","best", "Interpreter", "latex", "FontSize",16)