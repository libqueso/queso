clear all;
close all;
cd outputData
rawChain_ml;

fprintf(1,' Plotting KDE - theta_i, i=1..15 <press any key>\n');
for j=1:15
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,j),'function','pdf');
plot(xi,f,'-b','linewidth',3);
clear xi; clear f;
title(['\theta_{',num2str(j),'}: Kernel Density Estimation (raw chain)'],'fontsize',16);
set(gca,'FontSize',16);
ylabel('KDE','fontsize',16);
grid minor;
print(gcf,'-dpng',strcat('hysteretic_kde_theta',num2str(j)))
end


fprintf(1,' Plotting CDF - theta_i, i=1..15 <press any key>\n');
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,1),'function','cdf');
plot(xi,f,'-r','linewidth',3); hold on; clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,2),'function','cdf');
plot(xi,f,'-g','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,3),'function','cdf');
plot(xi,f,'-b','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,4),'function','cdf');
plot(xi,f,'-c','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,5),'function','cdf');
plot(xi,f,'-m','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,6),'function','cdf');
plot(xi,f,'-k','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,7),'function','cdf');
plot(xi,f,'--r','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,8),'function','cdf');
plot(xi,f,'--g','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,9),'function','cdf');
plot(xi,f,'--b','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,10),'function','cdf');
plot(xi,f,'--c','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,11),'function','cdf');
plot(xi,f,'--m','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,12),'function','cdf');
plot(xi,f,'--k','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,13),'function','cdf');
plot(xi,f,'-.r','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,14),'function','cdf');
plot(xi,f,'-.g','linewidth',3); clear xi; clear f;
[f,xi] = ksdensity(ip_ml_13_rawChain_unified(:,15),'function','cdf');
plot(xi,f,'-.b','linewidth',3); clear xi; clear f;
title('\theta_i: Cumulative Distribution Function (raw chain)','fontsize',16);
h=legend('\theta_1','\theta_2','\theta_3','\theta_4','\theta_5',...
         '\theta_6','\theta_7','\theta_8','\theta_9','\theta_{10}',...
         '\theta_{11}','\theta_{12}','\theta_{13}','\theta_{14}',...
         '\theta_{15}','location','EastOutside');
grid minor;
hold off; 
axis([-1.5 1.5 0 1]);
print -dpng hysteretic_cdf_thetas


fprintf(1,' Plotting autocorrelation - theta_i, i=1..15 <press any key>\n');
nlags=20;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,1), nlags,0);
plot(ACF1, lags1,'-r','linewidth',3); hold on; clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,2), nlags,0);
plot(ACF1, lags1,'-g','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,3), nlags,0);
plot(ACF1, lags1,'-b','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,4), nlags,0);
plot(ACF1, lags1,'-c','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,5), nlags,0);
plot(ACF1, lags1,'-m','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,6), nlags,0);
plot(ACF1, lags1,'-k','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,7), nlags,0);
plot(ACF1, lags1,'--r','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,8), nlags,0);
plot(ACF1, lags1,'--g','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,9), nlags,0);
plot(ACF1, lags1,'--b','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,10), nlags,0);
plot(ACF1, lags1,'--c','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,11), nlags,0);
plot(ACF1, lags1,'--m','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,12), nlags,0);
plot(ACF1, lags1,'--k','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,13), nlags,0);
plot(ACF1, lags1,'-.r','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,14), nlags,0);
plot(ACF1, lags1,'-.g','linewidth',3); clear ACF1; clear lags1;
[ACF1,lags1,bounds] = autocorr(ip_ml_13_rawChain_unified(:,15), nlags,0);
plot(ACF1, lags1,'-.b','linewidth',3); clear ACF1; clear lags1;
title('\theta_i: Autocorrelation (raw chain)','fontsize',16);
h=legend('\theta_1','\theta_2','\theta_3','\theta_4','\theta_5',...
         '\theta_6','\theta_7','\theta_8','\theta_9','\theta_{10}',...
         '\theta_{11}','\theta_{12}','\theta_{13}','\theta_{14}',...
         '\theta_{15}','location','EastOutside');
grid minor;
hold off; 
axis([0.3 1 0 20]);
print -dpng hysteretic_autocorr_thetas
clf;

% scatter plots - parameters
% plotmatrix(ip_ml_13_rawChain_unified)
gplotmatrix(ip_ml_13_rawChain_unified)
print -dpng hysteretic_scatter_thetas


% Likelihood plots
fprintf(1,' Plotting KDE - log(likelihood) and likelihood <press any key>\n');
[f,xi] = ksdensity(exp(ip_ml_13_rawLogLikelihood_unified),'function','pdf');
plot(xi,f,'-k','linewidth',3); clear xi; clear f;
title('KDE of Likelihood (f({\bfy}|{\bf\theta})) - last level','fontsize',16);
grid minor;
set(gca,'FontSize',16);
print -dpng hysteretic_kde_likelihood

[f,xi] = ksdensity(ip_ml_13_rawLogLikelihood_unified,'function','pdf');
plot(xi,f,'-k','linewidth',3); clear xi; clear f;
title('KDE of the Logarithm of the Likelihood (log(f({\bfy}|{\bf\theta}))) - last level','fontsize',16);
grid minor;
set(gca,'FontSize',16);
print -dpng hysteretic_kde_log_likelihood


cd ..
