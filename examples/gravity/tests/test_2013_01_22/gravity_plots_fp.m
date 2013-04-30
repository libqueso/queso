clear all;
cd outputData/
sfp_gravity_qoi_seq
sfp_gravity_p_seq;


% Chain plot --------------------------------------------------------------
fprintf(1,' Plotting chain  <press any key>\n');
plot(fp_mc_QoiSeq_unified,'m');
ylabel('QoI','fontname', 'Times', 'fontsize',20);
xlabel('Number of positions','fontname', 'Times', 'fontsize',20);
title('MC Chain Positions','fontname', 'Times', 'fontsize',20);
%a=axis;
%axis([a(1) a(2) .024 .028]);
set(gca,'FontSize',16);
print -dpng sfp_gravity_chain_pos.png
pause;
clf;

% Histogram plots ---------------------------------------------------------
fprintf(1,' Plotting histogram  <press any key>\n');
nbins=100;
hist(fp_mc_QoiSeq_unified,nbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','m','EdgeColor','k');%
title('QoI Histogram (nbins=100)','fontname', 'Times', 'fontsize',20);
xlabel('Distance traveled (m)','fontname', 'Times', 'fontsize',20);
ylabel('Frequency','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng sfp_gravity_hist.png
pause;
clf;


% KDE plots ---------------------------------------------------------------
fprintf(1,' Plotting KDE <press any key>\n');
[f,xi] = ksdensity(fp_mc_QoiSeq_unified,'function','pdf');
plot(xi,f,'-m','linewidth',3)
title('QoI Kernel Density Estimation ','fontname', 'Times', 'fontsize',20);
xlabel('Distance traveled (m)','fontname', 'Times', 'fontsize',20);
ylabel('KDE','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng sfp_gravity_kde.png
pause;
clf;

% CDF plots ---------------------------------------------------------------
fprintf(1,' Plotting CDF <press any key>\n');
[f,xi] = ksdensity(fp_mc_QoiSeq_unified,'function','cdf');
plot(xi,f,'-m','linewidth',3)
title('QoI Cumulative Distribution Function ','fontname', 'Times', 'fontsize',20);
xlabel('Distance traveled (m)','fontname', 'Times', 'fontsize',20);
ylabel('CDF','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng sfp_gravity_cdf.png
pause;
clf;


% Autocorrelation plots ---------------------------------------------------
%Let's use the autocorr function provided in
% http://sfb649.wiwi.hu-berlin.de/quantnet/index.php?p=show&id=900

fprintf(1,' Plotting autocorrelation  <press any key>\n');
cd ../
nlags=10;
[ACF, lags, bounds] = autocorr_local(fp_mc_QoiSeq_unified, nlags, 0);
plot(lags,ACF,'mo-','linewidth',3);
ylabel('Autocorrelation for QoI = d','fontname', 'Times', 'fontsize',20);
xlabel('Lag','fontname', 'Times', 'fontsize',20);
title('QoI Autocorrelation','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng outputData/sfp_gravity_autocorrelation.png
pause;
clf;
close;


% Autocorrelation plots ---------------------------------------------------
% If Matlab has Econometrics Toolbox, then you may use its own autocorr 
% function
% 
% fprintf(1,' Plotting autocorrelation  <press any key>\n');
% nlags=10;
% [ACF, lags, bounds] = autocorr(fp_mc_QoiSeq_unified, nlags, 0);
% plot(lags,ACF,'mo-','linewidth',3);
% ylabel('Autocorrelation for QoI = d','fontname', 'Times', 'fontsize',20);
% xlabel('Lag','fontname', 'Times', 'fontsize',20);
% title('QoI Autocorrelation','fontname', 'Times', 'fontsize',20);
% grid on;
% set(gca,'FontSize',16);
% print -dpng sfp_gravity_autocorrelation.png
% pause;
% clf;



% Covariance matrix -------------------------------------------------------
fprintf(1, '\n Covariance matrix:  <press any key>\n');
X=[fp_mc_ParamSeq_unified fp_mc_QoiSeq_unified];
cov_d = cov(X);
fprintf(1,' cov(p,QoI) =\n\t [ %3.3e \t %3.3e ] \n\t [ %3.3e \t %3.3e ]\n',cov_d(1,1), cov_d(1,2), cov_d(2,1), cov_d(2,2));


% Correlation matrix ------------------------------------------------------
fprintf(1, '\n Correlation matrix:  <press any key>\n');
X=[fp_mc_ParamSeq_unified fp_mc_QoiSeq_unified];
corr_d = corrcoef(X);
fprintf(1,' corr(p,QoI) =\n\t [ %3.3e \t %3.3e ] \n\t [ %3.3e \t %3.3e ]\n',corr_d(1,1), corr_d(1,2), corr_d(2,1), corr_d(2,2));

% -------------------------------------------------------------------------

