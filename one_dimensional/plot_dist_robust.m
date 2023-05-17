%% Load Data
clear;clc;
d = [5,10,15,20];
X ={};
for i=1:4
    matching_string = sprintf("*_d%i_*",d(i));
    myFiles_norm = dir(fullfile("Normal_1d",matching_string));
    myFiles_exp = dir(fullfile("Exponential_1d",matching_string));
    for k = 1:length(myFiles_norm)
        filename_norm = strcat(myFiles_norm(k).folder,'/',myFiles_norm(k).name);
        filename_exp = strcat(myFiles_exp(k).folder,'/',myFiles_exp(k).name);
        X{i,k} = load(filename_norm);
        Y{i,k} = load(filename_exp);
    end
end

%% Compute Means
[n_pts,n_mc]  = size(X{1,1}.dist_v_n);
n_space = X{1,1}.n_space;
for i =1:4
    dist_v_n_norm = zeros(n_pts,n_mc);
    dist_v_n_exp = zeros(n_pts,n_mc);
    for k=1:10
        dist_v_n_norm = dist_v_n_norm + X{i,k}.dist_v_n;
        dist_v_n_exp = dist_v_n_exp + Y{i,k}.dist_v_n;
    end
    norm_mean(:,i) = mean(dist_v_n_norm/10,2, "omitnan");
    exp_mean(:,i) = mean(dist_v_n_exp/10,2, "omitnan");
end
%% Compute std
for i=1:4
    dist_v_n_norm = zeros(n_pts,n_mc);
    dist_v_n_exp = zeros(n_pts,n_mc);
    for k=1:10
        dist_v_n_norm = dist_v_n_norm + (X{i,k}.dist_v_n - norm_mean(:,i)).^2;
        dist_v_n_exp = dist_v_n_exp + (Y{i,k}.dist_v_n - exp_mean(:,i)).^2;
    end
    norm_var(:,i) = sum(dist_v_n_norm/10,2, "omitnan")/(10*n_mc - 1);
    exp_var(:,i) = sum(dist_v_n_exp/10,2, "omitnan")/(10*n_mc - 1);
end
normstd = norm_var.^0.5;
expstd = exp_var.^0.5;
normstderr = normstd/(n_mc*10)^0.5;
expstderr = expstd/(n_mc*10)^0.5;
figure;
a = axes;
errorbar(n_space, norm_mean(:,1), normstd(:,1), "Linewidth", 2, "Color", "#0072BD");
hold on
errorbar(n_space, exp_mean(:,1), expstd(:,1), "Linewidth", 2, "Color", "#D95319");
errorbar(n_space, norm_mean(:,4), normstd(:,4),"-." ,"Linewidth", 2, "Color", "#0072BD");
hold on
errorbar(n_space, exp_mean(:,4), expstd(:,4), "-.","Linewidth", 2, "Color", "#D95319");
ylim([2e-3 2e-1])
a.XScale = "log";
a.YScale = "log";
xlabel("Number of Samples $n$", "Interpreter","latex", "FontSize",16)
ylabel("$d_{TV}(\hat{p}_(y|\mathbf{x}), p_(y|\mathbf{x}))$", "Interpreter","latex", "FontSize",16)
legend("$\mathbf{x} \sim\mathcal{N}(0,I_5)$","$\mathbf{x}\sim \bigotimes_5 \mathrm{Lap}\left( 0,1\right)$", "$\mathbf{x} \sim\mathcal{N}(0,I_{20})$","$\mathbf{x}\sim\bigotimes_{20}\mathrm{Lap}\left( 0,1\right), k =20$","Interpreter", "Latex","FontSize",16)