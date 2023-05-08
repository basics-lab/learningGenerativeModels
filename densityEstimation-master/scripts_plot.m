clear all;
load('fig1_subfig1.mat')
err_b_mean = mean(err_b, 2);
err_b_std = std(err_b, 0, 2);
err_sigma_mean = mean(err_sigma, 2);
err_sigma_std = std(err_sigma, 0, 2);
sqrt_kl_mean = mean(sqrt_kl, 2); 
sqrt_kl_std = std(sqrt_kl, 0, 2);
subplot(1,3,1);
plot(num_samples/1e5, err_b_mean, '^-', 'linewidth', 1.0, 'color', [0, 0.44, 0.74]);
hold on
plot(num_samples/1e5, err_sigma_mean, 'o-', 'linewidth', 1.0, 'color', [0.847, 0.32, 0.094]);
plot(num_samples/1e5, sqrt_kl_mean, '*-', 'linewidth', 1.0, 'color', [0.52, 0.76, 0.21]);
xlabel('Number of Samples (x10^5)', 'FontSize',9);
ylabel('Error', 'FontSize', 9);
set(gca,'FontSize',9)
clear all;
load('fig1_subfig2.mat')
err_b_mean = mean(err_b, 2);
err_b_std = std(err_b, 0, 2);
err_sigma_mean = mean(err_sigma, 2);
err_sigma_std = std(err_sigma, 0, 2);
sqrt_kl_mean = mean(sqrt_kl, 2); 
sqrt_kl_std = std(sqrt_kl, 0, 2);
subplot(1,3,2);
plot(dims, err_b_mean, '^-', 'linewidth', 1.0, 'color', [0, 0.44, 0.74]);
hold on
plot(dims, err_sigma_mean, 'o-', 'linewidth', 1.0, 'color', [0.847, 0.32, 0.094]);
plot(dims, sqrt_kl_mean, '*-', 'linewidth', 1.0, 'color', [0.52, 0.76, 0.21]);
xlabel('Dimension ({\itd})', 'FontSize',9);
ylabel('Error', 'FontSize',9);
set(gca,'FontSize',9)
axis([5,25,0,0.25]);
clear all;
load('fig1_subfig3.mat')
err_b_mean = mean(err_b, 2);
err_b_std = std(err_b, 0, 2);
err_sigma_mean = mean(err_sigma, 2);
err_sigma_std = std(err_sigma, 0, 2);
sqrt_kl_mean = mean(sqrt_kl, 2); 
sqrt_kl_std = std(sqrt_kl, 0, 2);
subplot(1,3,3);
h1 = plot(kappas, err_b_mean, '^-', 'linewidth', 1.0, 'color', [0, 0.44, 0.74]);
hold on
h2 = plot(kappas, err_sigma_mean, 'o-', 'linewidth', 1.0, 'color', [0.847, 0.32, 0.094]);
h3 = plot(kappas, sqrt_kl_mean, '*-', 'linewidth', 1.0, 'color', [0.52, 0.76, 0.21]);
xlabel('Condition Number (\kappa)', 'FontSize',9);
ylabel('Error', 'FontSize',9);
set(gca,'FontSize',9)
leg = legend([h1, h2, h3], {'$||\widehat{b}-b^*||_2/||W^*||_F$', 
    '$||\widehat{\Sigma}-W^*W^{*T}||_2/||W^*||_F^2$',
    '$\sqrt{KL(N(\widehat{b}, \widehat{\Sigma})||N(b^*, W^*W^{*T}))/2}$'});
set(leg,'Interpreter','latex', 'FontSize',12, 'Orientation', 'horizontal');
axis([1,25,0,0.12]);