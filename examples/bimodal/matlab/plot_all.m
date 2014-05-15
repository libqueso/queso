clear all;
%cd ../test/test_2013_12_02/outputData
cd outputData
rawChain_ml


% KDE estimation
fprintf(1,' Plotting KDE - theta <press any key>\n');
% [f,xi] = ksdensity(ip_ml_0_rawChain_unified(:,1),'function','pdf');
% plot(xi,f,'-m','linewidth',3);
% hold on;  clear xi; clear f;
[f,xi] = ksdensity(ip_ml_1_rawChain_unified(:,1),'function','pdf');
plot(xi,f,'-r','linewidth',3);
hold on;  clear xi; clear f;
[f,xi] = ksdensity(ip_ml_2_rawChain_unified(:,1),'function','pdf');
plot(xi,f,'-g','linewidth',3)
[f,xi] = ksdensity(ip_ml_3_rawChain_unified(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',3)
[f,xi] = ksdensity(ip_ml_4_rawChain_unified(:,1),'function','pdf');
plot(xi,f,'-c','linewidth',3)
% h=legend('Level0','Level 1','Level 2','Level 3','Level 4','location','northeast');
h=legend('Level 1','Level 2','Level 3','Level 4','location','northeast');
title('\theta Kernel Density Estimation (raw chain)','fontname', 'Times', 'fontsize',20);
set(gca,'FontSize',16);
grid minor;
a=axis;
axis([-50 200 a(3) a(4)]);
hold off;
print -dpng bimodal_kde_rawchain
pause;
clf; 

% CDF estimation
fprintf(1,' Plotting CDF - theta <press any key>\n');
[f,xi] = ksdensity(ip_ml_1_rawChain_unified(:,1),'function','cdf');
plot(xi,f,'-r','linewidth',3);
hold on;  clear xi; clear f;
[f,xi] = ksdensity(ip_ml_2_rawChain_unified(:,1),'function','cdf');
plot(xi,f,'-g','linewidth',3)
clear xi; clear f;
[f,xi] = ksdensity(ip_ml_3_rawChain_unified(:,1),'function','cdf');
plot(xi,f,'-b','linewidth',3)
clear xi; clear f;
[f,xi] = ksdensity(ip_ml_4_rawChain_unified(:,1),'function','cdf');
plot(xi,f,'-c','linewidth',3)
clear xi; clear f;
h=legend('Level 1','Level 2','Level 3','Level 4','location','southeast');
title('\theta Cumulative Distribution Function (raw chain)','fontname', 'Times', 'fontsize',20);set(gca,'FontSize',16);
grid minor;
a=axis;
axis([-50 200 a(3) a(4)]);
hold off;
print -dpng bimodal_cdf_rawchain
pause;

% AUTOCORRELATION PLOTS
fprintf(1,' Plotting autocorrelation - theta <press any key>\n');
nlags=20;
i=1;
[ACF0,lags0,bounds]= autocorr(ip_ml_0_rawChain_unified(:,i), nlags, 0);
[ACF1,lags1,bounds]= autocorr(ip_ml_1_rawChain_unified(:,i), nlags, 0);
[ACF2,lags2,bounds]= autocorr(ip_ml_2_rawChain_unified(:,i), nlags, 0);
[ACF3,lags3,bounds]= autocorr(ip_ml_3_rawChain_unified(:,i), nlags, 0);
[ACF4,lags4,bounds]= autocorr(ip_ml_4_rawChain_unified(:,i), nlags, 0);
plot(ACF1,lags1,'-r',ACF2,lags2,'-g',ACF3,lags3,'-b',ACF4,lags4,'-c','linewidth',3);
h=legend('Level 1','Level 2','Level 3','Level 4','location','northeast');
set(gca,'FontSize',16);
grid minor;
ylabel('Autocorrelation', 'fontsize',16);
xlabel('Lag','fontname', 'Times', 'fontsize',16);
title('Autocorrelation - Parameter \theta','fontname', 'Times','fontsize',20);
print -dpng bimodal_autocorrelation_rawchain
pause;
clf; 

%cd ../../../matlab
cd ..
