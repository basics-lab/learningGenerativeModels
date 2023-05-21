clear;clc;
myFiles = dir(fullfile("data",'*.mat'));
X ={};
n_runs = length(myFiles);
for k = 1:length(myFiles)
    filename = strcat(myFiles(k).folder,'/',myFiles(k).name);
    X{k} = load(filename);
end
N = X{1}.N;
n_mc = X{1}.n_mc;
b_space = linspace(-5,5,N);
dist_v_b = zeros(N,8);
dist_v_b_wu = zeros(N,8);
dist_v_b_wu2 = zeros(N,8);

for i=1:n_runs
    dist_v_b = dist_v_b + X{i}.dist_v_b;
    dist_v_b_wu = dist_v_b_wu + X{i}.dist_v_b_wu;
    dist_v_b_wu2 = dist_v_b_wu2 + X{i}.dist_v_b_wu2;
end
mean_dist = dist_v_b/n_runs;
mean_dist_wu = dist_v_b_wu/n_runs;
mean_dist_wu2 = dist_v_b_wu2/n_runs;
f = figure;
t = tiledlayout(1,1);
a1 = axes(t);
p1 = plot(a1,b_space, mean_dist(:,1), "LineWidth", 3);
hold on;
p2 = plot(a1,b_space(5:end), mean_dist_wu(5:end,1), "LineWidth", 3, "Color", "#D95319");
%plot(b_space(11:end), mean_dist_wu2(11:end,1), "Linewidth", 2);
a1.YScale = 'log';
xlim([-4 5])
ylim([1e-4 0.3])
xlabel("Bias Value $b$", "Interpreter","latex", "FontSize",16)
ylabel("$d\left((\hat{w},\hat{\Sigma}), ( w^*,\Sigma^*) \right)$", "Interpreter","latex", "FontSize",16)
a2 = axes(t);
hold on
l2_wu = mean_dist_wu(:,4).^2 + mean_dist_wu(:,3).^2; %mean_dist_wu(:,7).^2 + ;
l2 =  mean_dist(:,4).^2 + mean_dist(:,3).^2; %mean_dist(:,7).^2;
p3 = plot(a2,b_space(4:end), l2(4:end), "--", "LineWidth", 1);
p4 = plot(a2,b_space(7:end), l2_wu(7:end), "--", "LineWidth", 1 );
xlim([-4 5])
ylim([2e-4 0.3])
a2.YScale = 'log';
a2.XAxis.Visible = 'off';
a2.YAxisLocation = 'right';
a2.Color = 'none';
a1.Box = 'on';
a2.Box = 'off';
ylabel("$||\hat{\Sigma} - \Sigma^*||^2_F + ||\hat{w} - w^*||^2$", "Interpreter","latex", "FontSize",16)
legend([p1,p2,p3,p4], ... 
    "Our TV", ...
    "Wu et al. TV",...
    "Our parameter MSE",...
    "Wu et al. parameter MSE", ...
    "Location","best", "Interpreter", "latex", "FontSize",14, 'NumColumns',1);