clear all;
cd ../tests/test_2013_11_22/outputData_1mode
rawChain_ml

% Level 0
fprintf(1,'Scatter plots and histograms of raw chains - Level 0 <press any key>\n');
plotmatrix(ip_ml_0_rawChain_unified, '+b')
set(gca,'fontsize',20); 
xlabel('\theta_1                       \theta_2                        \theta_3','fontsize',16);
ylabel('\theta_3                       \theta_2                        \theta_1','fontsize',16);
title('Scatter plots and histograms, Level 0 - 1 mode')
print -dpng modal_1_mode_level_0
pause;
clf; 

% Level 1
fprintf(1,'Scatter plots and histograms of raw chains - Level 1 <press any key>\n');
plotmatrix(ip_ml_1_rawChain_unified, '+b')
set(gca,'fontsize',20); 
xlabel('\theta_1                       \theta_2                        \theta_3','fontsize',16);
ylabel('\theta_3                       \theta_2                        \theta_1','fontsize',16);
title('Scatter plots and histograms, Level 1 - 1 mode')
print -dpng modal_1_mode_level_1
pause;
clf; 

% Level 2
fprintf(1,' Scatter plots and histograms of raw chains - Level 2 <press any key>\n');
plotmatrix(ip_ml_2_rawChain_unified, '+b')
set(gca,'fontsize',20); 
xlabel('\theta_1                       \theta_2                        \theta_3','fontsize',16);
ylabel('\theta_3                       \theta_2                        \theta_1','fontsize',16);
title('Scatter plots and histograms, Level 2 - 1 mode')
print -dpng modal_1_mode_level_2
pause;
clf; 

% Level 3
fprintf(1,' Scatter plots and histograms of raw chains - Level 3 <press any key>\n');
plotmatrix(ip_ml_3_rawChain_unified, '+b')
set(gca,'fontsize',20); 
xlabel('\theta_1                       \theta_2                        \theta_3','fontsize',16);
ylabel('\theta_3                       \theta_2                        \theta_1','fontsize',16);
title('Scatter plots and histograms, Level 3 - 1 mode')
print -dpng modal_1_mode_level_3
pause;
clf; 

% Level 4
fprintf(1,' Scatter plots and histograms of raw chains - Level 4 <press any key>\n');
plotmatrix(ip_ml_4_rawChain_unified, '+b')
set(gca,'fontsize',20); 
xlabel('\theta_1                       \theta_2                        \theta_3','fontsize',16);
ylabel('\theta_3                       \theta_2                        \theta_1','fontsize',16);
title('Scatter plots and histograms, Level 4 - 1 mode')
print -dpng modal_1_mode_level_4
pause;
clf; 

% Level 5
fprintf(1,' Scatter plots and histograms of raw chains - Level 5 <press any key>\n');
plotmatrix(ip_ml_5_rawChain_unified, '+b')
set(gca,'fontsize',20); 
xlabel('\theta_1                       \theta_2                        \theta_3','fontsize',16);
ylabel('\theta_3                       \theta_2                        \theta_1','fontsize',16);
title('Scatter plots and histograms, Level 5 - 1 mode')
print -dpng modal_1_mode_level_5
pause;
clf; 

% Level 6
fprintf(1,' Scatter plots and histograms of raw chains - Level 6 <press any key>\n');
plotmatrix(ip_ml_6_rawChain_unified, '+b')
set(gca,'fontsize',20); 
xlabel('\theta_1                       \theta_2                        \theta_3','fontsize',16);
ylabel('\theta_3                       \theta_2                        \theta_1','fontsize',16);
title('Scatter plots and histograms, Level 6 - 1 mode')
print -dpng modal_1_mode_level_6
pause;
clf; 

% Level 7
fprintf(1,' Scatter plots and histograms of raw chains - Level 7 <press any key>\n');
plotmatrix(ip_ml_7_rawChain_unified, '+b')
set(gca,'fontsize',20); 
xlabel('\theta_1                       \theta_2                        \theta_3','fontsize',16);
ylabel('\theta_3                       \theta_2                        \theta_1','fontsize',16);
title('Scatter plots and histograms, Level 7 - 1 mode')
print -dpng modal_1_mode_level_7


% KDE estimation
% theta 1
fprintf(1,' Plotting KDE - theta_1 <press any key>\n');
[f,xi] = ksdensity(ip_ml_1_rawChain_unified(:,1),'function','pdf');
plot(xi,f,'-r','linewidth',3);
hold on;  clear xi; clear f;
[f,xi] = ksdensity(ip_ml_2_rawChain_unified(:,1),'function','pdf');
plot(xi,f,'-g','linewidth',3)
[f,xi] = ksdensity(ip_ml_3_rawChain_unified(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',3)
[f,xi] = ksdensity(ip_ml_4_rawChain_unified(:,1),'function','pdf');
plot(xi,f,'-c','linewidth',3)
[f,xi] = ksdensity(ip_ml_5_rawChain_unified(:,1),'function','pdf');
plot(xi,f,'-m','linewidth',3)
[f,xi] = ksdensity(ip_ml_6_rawChain_unified(:,1),'function','pdf');
plot(xi,f,'-k','linewidth',3)
[f,xi] = ksdensity(ip_ml_7_rawChain_unified(:,1),'function','pdf');
plot(xi,f,'--r','linewidth',3)
h=legend('Level 1','Level 2','Level 3','Level 4','Level 5','Level 6','Level 7','location','northeast');
title('\theta_1 Kernel Density Estimation (raw chain)','fontname', 'Times', 'fontsize',20);
set(gca,'FontSize',16);
grid on;
axis([-1 5 0 .6]);
hold off;
print -dpng modal_1_mode_kde_theta1
pause;
clf; 
 
%theta2
fprintf(1,' Plotting KDE - theta_2 <press any key>\n'); 
[f,xi] = ksdensity(ip_ml_1_rawChain_unified(:,2),'function','pdf');
plot(xi,f,'-r','linewidth',3);
hold on;  clear xi; clear f;
[f,xi] = ksdensity(ip_ml_2_rawChain_unified(:,2),'function','pdf');
plot(xi,f,'-g','linewidth',3)
[f,xi] = ksdensity(ip_ml_3_rawChain_unified(:,2),'function','pdf');
plot(xi,f,'-b','linewidth',3)
[f,xi] = ksdensity(ip_ml_4_rawChain_unified(:,2),'function','pdf');
plot(xi,f,'-c','linewidth',3)
[f,xi] = ksdensity(ip_ml_5_rawChain_unified(:,2),'function','pdf');
plot(xi,f,'-m','linewidth',3)
[f,xi] = ksdensity(ip_ml_6_rawChain_unified(:,2),'function','pdf');
plot(xi,f,'-k','linewidth',3)
[f,xi] = ksdensity(ip_ml_7_rawChain_unified(:,2),'function','pdf');
plot(xi,f,'--r','linewidth',3)
h=legend('Level 1','Level 2','Level 3','Level 4','Level 5','Level 6','Level 7','location','northeast');
title('\theta_2 Kernel Density Estimation (raw chain)','fontname', 'Times', 'fontsize',20);
set(gca,'FontSize',16);grid on;
axis([-1 4 0 1]);
hold off; 
print -dpng modal_1_mode_kde_theta2
pause;
clf; 



%theta 3 
fprintf(1,' Plotting KDE - theta_3 <press any key>\n');
[f,xi] = ksdensity(ip_ml_1_rawChain_unified(:,3),'function','pdf');
plot(xi,f,'-r','linewidth',3);
hold on;  clear xi; clear f;
[f,xi] = ksdensity(ip_ml_2_rawChain_unified(:,3),'function','pdf');
plot(xi,f,'-g','linewidth',3)
[f,xi] = ksdensity(ip_ml_3_rawChain_unified(:,3),'function','pdf');
plot(xi,f,'-b','linewidth',3)
[f,xi] = ksdensity(ip_ml_4_rawChain_unified(:,3),'function','pdf');
plot(xi,f,'-c','linewidth',3)
[f,xi] = ksdensity(ip_ml_5_rawChain_unified(:,3),'function','pdf');
plot(xi,f,'-m','linewidth',3)
[f,xi] = ksdensity(ip_ml_6_rawChain_unified(:,3),'function','pdf');
plot(xi,f,'-k','linewidth',3)
[f,xi] = ksdensity(ip_ml_7_rawChain_unified(:,3),'function','pdf');
plot(xi,f,'--r','linewidth',3)
h=legend('Level 1','Level 2','Level 3','Level 4','Level 5','Level 6','Level 7','location','northeast');
title('\theta_3 Kernel Density Estimation (raw chain)','fontname', 'Times', 'fontsize',20);
set(gca,'FontSize',16);grid on;
axis([-0.2 .8 0 30]);
hold off; 
print -dpng modal_1_mode_kde_theta3
pause;
clf; 

% Target pdf
fprintf(1,' Plotting KDE - target PDF <press any key>\n');
%[f,xi] = ksdensity(ip_ml_0_rawLogTarget_unified,'function','pdf');
%plot(xi,f,'-k','linewidth',3)
%hold on;  clear xi; clear f;
[f,xi] = ksdensity(ip_ml_1_rawLogTarget_unified,'function','pdf');
plot(xi,f,'-r','linewidth',3);
hold on;  clear xi; clear f;
[f,xi] = ksdensity(ip_ml_2_rawLogTarget_unified,'function','pdf');
plot(xi,f,'-g','linewidth',3)
[f,xi] = ksdensity(ip_ml_3_rawLogTarget_unified,'function','pdf');
plot(xi,f,'-b','linewidth',3)
[f,xi] = ksdensity(ip_ml_4_rawLogTarget_unified,'function','pdf');
plot(xi,f,'-c','linewidth',3)
[f,xi] = ksdensity(ip_ml_5_rawLogTarget_unified,'function','pdf');
plot(xi,f,'-m','linewidth',3)
[f,xi] = ksdensity(ip_ml_6_rawLogTarget_unified,'function','pdf');
plot(xi,f,'--k','linewidth',3)
[f,xi] = ksdensity(ip_ml_7_rawLogTarget_unified,'function','pdf');
plot(xi,f,'--r','linewidth',3)
h=legend('Level 1','Level 2','Level 3','Level 4','Level 5','Level 6','Level 7','location','northwest');
title('Target PDF Kernel Density Estimation','fontname', 'Times', 'fontsize',20);
hold off; grid on;set(gca,'FontSize',16);
axis([0 15 0 0.6]);
print -dpng modal_1_mode_kde_target
pause;
clf; 

% AUTOCORRELATION PLOTS
nlags=20;
%Theta 1
i=1;
[ACF0,lags,bounds]= autocorr(ip_ml_0_rawChain_unified(:,i), nlags, 0);
[ACF1,lags,bounds]= autocorr(ip_ml_1_rawChain_unified(:,i), nlags, 0);
[ACF2,lags,bounds]= autocorr(ip_ml_2_rawChain_unified(:,i), nlags, 0);
[ACF3,lags,bounds]= autocorr(ip_ml_3_rawChain_unified(:,i), nlags, 0);
[ACF4,lags,bounds]= autocorr(ip_ml_4_rawChain_unified(:,i), nlags, 0);
[ACF5,lags,bounds]= autocorr(ip_ml_5_rawChain_unified(:,i), nlags, 0);
[ACF6,lags,bounds]= autocorr(ip_ml_6_rawChain_unified(:,i), nlags, 0);
[ACF7,lags,bounds]= autocorr(ip_ml_7_rawChain_unified(:,i), nlags, 0);
plot(ACF0,lags,'-k', ACF1,lags,'-r',ACF2,lags,'-g',ACF3,lags,'-b',ACF4,lags,'-c',ACF5,lags,'-m',ACF6,lags,'--k',ACF7,lags,'--r','linewidth',3);
h=legend('Level 1','Level 2','Level 3','Level 4','Level 5','Level 6','Level 7','location','northeast');
set(gca,'FontSize',16);
grid minor;
title('Autocorrelation - Parameter \theta_1','fontname', 'Times', 'fontsize',20);
print -dpng modal_1_mode_autocorrelation_theta1
pause;
clf; 

% AUTOCORRELATION PLOTS
%Theta 2
i=2;
[ACF0,lags,bounds]= autocorr(ip_ml_0_rawChain_unified(:,i), nlags, 0);
[ACF1,lags,bounds]= autocorr(ip_ml_1_rawChain_unified(:,i), nlags, 0);
[ACF2,lags,bounds]= autocorr(ip_ml_2_rawChain_unified(:,i), nlags, 0);
[ACF3,lags,bounds]= autocorr(ip_ml_3_rawChain_unified(:,i), nlags, 0);
[ACF4,lags,bounds]= autocorr(ip_ml_4_rawChain_unified(:,i), nlags, 0);
[ACF5,lags,bounds]= autocorr(ip_ml_5_rawChain_unified(:,i), nlags, 0);
[ACF6,lags,bounds]= autocorr(ip_ml_6_rawChain_unified(:,i), nlags, 0);
[ACF7,lags,bounds]= autocorr(ip_ml_7_rawChain_unified(:,i), nlags, 0);
plot(ACF0,lags,'-k', ACF1,lags,'-r',ACF2,lags,'-g',ACF3,lags,'-b',ACF4,lags,'-c',ACF5,lags,'-m',ACF6,lags,'--k',ACF7,lags,'--r','linewidth',3);
h=legend('Level 0','Level 1','Level 2','Level 3','Level 4','Level 5','Level 6','Level 7','location','northeast');
set(gca,'FontSize',16);
grid minor;
title('Autocorrelation - Parameter \theta_2','fontname', 'Times', 'fontsize',20);
print -dpng modal_1_mode_autocorrelation_theta2
pause;
clf; 


%Theta 3
i=3;
[ACF0,lags,bounds]= autocorr(ip_ml_0_rawChain_unified(:,i), nlags, 0);
[ACF1,lags,bounds]= autocorr(ip_ml_1_rawChain_unified(:,i), nlags, 0);
[ACF2,lags,bounds]= autocorr(ip_ml_2_rawChain_unified(:,i), nlags, 0);
[ACF3,lags,bounds]= autocorr(ip_ml_3_rawChain_unified(:,i), nlags, 0);
[ACF4,lags,bounds]= autocorr(ip_ml_4_rawChain_unified(:,i), nlags, 0);
[ACF5,lags,bounds]= autocorr(ip_ml_5_rawChain_unified(:,i), nlags, 0);
[ACF6,lags,bounds]= autocorr(ip_ml_6_rawChain_unified(:,i), nlags, 0);
[ACF7,lags,bounds]= autocorr(ip_ml_7_rawChain_unified(:,i), nlags, 0);
plot(ACF0,lags,'-k', ACF1,lags,'-r',ACF2,lags,'-g',ACF3,lags,'-b',ACF4,lags,'-c',ACF5,lags,'-m',ACF6,lags,'--k',ACF7,lags,'--r','linewidth',3);
h=legend('Level 0','Level 1','Level 2','Level 3','Level 4','Level 5','Level 6','Level 7','location','northeast');
set(gca,'FontSize',16);
grid minor;
title('Autocorrelation - Parameter \theta_3','fontname', 'Times', 'fontsize',20);
print -dpng modal_1_mode_autocorrelation_theta3
pause;
clf; 



cd ../matlab

