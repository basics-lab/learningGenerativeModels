%% Load files
clear;clc;
myFiles = dir(fullfile("data_may_14_attempt2",'*.mat'));
X ={};
n_runs = length(myFiles);
for k = 1:length(myFiles)
    filename = strcat(myFiles(k).folder,'/',myFiles(k).name);
    X{k} = load(filename);
end
N = X{1}.N;
n_mc = X{1}.n_mc;
kappa_space = X{1}.kappa;
dist_v_cond = zeros(6,N);
dist_v_cond_wu = zeros(6,N);
%% Mean
for i=1:n_runs
    dist_v_cond = dist_v_cond + squeeze(mean(X{i}.dist_v_cond, 1, "omitnan"));
    dist_v_cond_wu = dist_v_cond_wu + squeeze(mean(X{i}.dist_v_cond_wu, "omitnan"));
end
mean_dist = dist_v_cond/n_runs;
mean_dist_wu = dist_v_cond_wu/n_runs;
%% Variance
var_v_cond = zeros(6,N);
var_v_cond_wu = zeros(6,N);
for i=1:n_runs
    for j=1:N
        var_v_cond(:,j) = var_v_cond(:,j) + mean((squeeze(X{i}.dist_v_cond(:,:,j))' - mean_dist(:,j)).^2, 2,"omitnan");
        var_v_cond_wu(:,j) = var_v_cond_wu(:,j) + mean((squeeze(X{i}.dist_v_cond_wu(:,:,j))' - mean_dist_wu(:,j)).^2, 2,"omitnan");
    end
end
emp_std = (var_v_cond/n_runs).^0.5;
emp_std_wu = (var_v_cond_wu/n_runs).^0.5;
std_err = emp_std/sqrt(n_mc*n_runs);
std_err_wu = emp_std/sqrt(n_mc*n_runs);
figure;
r = [1:3:20, 20];
errorbar(kappa_space(r).^2, mean_dist(1,r), std_err(1,r), "LineWidth", 2);
hold on;
errorbar(kappa_space(r).^2, mean_dist_wu(1,r), std_err_wu(1,r), "LineWidth", 2);
% x = [kappa_space, fliplr(kappa_space)];
% curve1 = [mean_dist_wu + 0.1*var_v_cond_wu, fliplr(mean_dist_wu - 0.1*var_v_cond_wu)];
% fill(x,curve1(1,:),[0.55, 0.81, 0.89]);
% hold on;
% curve2 = [mean_dist + 0.1*var_v_cond, fliplr(mean_dist - 0.1*var_v_cond)];
% fill(x,curve2(1,:),[1, 0.4, 0.2]);
%plot(kappa_space, mean_dist(1,:));
%plot(kappa_space, mean_dist_wu(1,:));
xlim([1 24^2])
ylim([0.012 0.032])
xlabel("Condition Number of $\Sigma^*$", "Interpreter","latex", "FontSize",16)
ylabel("$d\left((\hat{w},\hat{\Sigma}), ( w^*,\Sigma^*) \right)$", "Interpreter","latex", "FontSize",16)
legend("Our TV", ...
    "Wu et al. TV",  "Location","best", "Interpreter", "latex", "FontSize",16)