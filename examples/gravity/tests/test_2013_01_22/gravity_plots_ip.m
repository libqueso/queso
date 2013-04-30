clear all;
cd outputData/
sip_gravity_raw_chain
sip_gravity_filtered_chain


% Chain plots -------------------------------------------------------------
% RAW
fprintf(1,' Plotting chains - raw  <press any key>\n');
plot(ip_mh_rawChain_unified)
ylabel('\theta=g','fontname', 'Times', 'fontsize',20);
xlabel('Number of positions','fontname', 'Times', 'fontsize',20);
title('DRAM MCMC Chain Positions (raw)','fontname', 'Times', 'fontsize',20);
set(gca,'FontSize',16);
print -dpng sip_gravity_chain_pos_raw.png
pause;
clf;

% FILTERED 
fprintf(1,' Plotting chains - filtered  <press any key>\n');
plot(ip_mh_filtChain_unified,'r')
ylabel('\theta=g','fontname', 'Times', 'fontsize',20);
xlabel('Number of positions','fontname', 'Times', 'fontsize',20);
title('DRAM MCMC Chain Positions (filtered)','fontname', 'Times', 'fontsize',20);
a=axis;
axis([a(1) a(2) 8 10.5]);
set(gca,'FontSize',16);
print -dpng sip_gravity_chain_pos_filt.png
pause;
clf;

% Histogram plots ---------------------------------------------------------
% RAW
fprintf(1,' Plotting histogram - raw  <press any key>\n');
nbins=100;
hist(ip_mh_rawChain_unified,nbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','w');%
title('Parameter Histogram (raw chain, nbins=100)','fontname', 'Times', 'fontsize',20);
xlabel('Gravity (m/s^2)','fontname', 'Times', 'fontsize',20);
ylabel('Frequency','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng sip_gravity_hist_raw.png
pause;
clf;

%FILTERED 
fprintf(1,' Plotting histogram - filtered  <press any key>\n');
hist(ip_mh_filtChain_unified,nbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w');%
title('Parameter Histogram (filtered chain, nbins=100)','fontname', 'Times', 'fontsize',20);
xlabel('Gravity (m/s^2)','fontname', 'Times', 'fontsize',20);
ylabel('Frequency','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng sip_gravity_hist_filt.png
pause;
clf;

% KDE plots ---------------------------------------------------------------
% RAW
fprintf(1,' Plotting KDE - raw  <press any key>\n');
[f,xi] = ksdensity(ip_mh_rawChain_unified,'function','pdf');
plot(xi,f,'-b','linewidth',3)
title('Parameter Kernel Density Estimation (raw chain)','fontname', 'Times', 'fontsize',20);
xlabel('Gravity (m/s^2)','fontname', 'Times', 'fontsize',20);
ylabel('KDE','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng sip_gravity_kde_raw.png
pause;
clf;

%FILTERED
fprintf(1,' Plotting KDE - filtered  <press any key>\n');
[f,xi] = ksdensity(ip_mh_filtChain_unified,'function','pdf');
plot(xi,f,'-r','linewidth',3)
title('Parameter Kernel Density Estimation (filtered chain)','fontname', 'Times', 'fontsize',20);
xlabel('Gravity (m/s^2)','fontname', 'Times', 'fontsize',20);
ylabel('KDE','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng sip_gravity_kde_filt.png
pause;
clf;

% CDF plots ---------------------------------------------------------------
% RAW
fprintf(1,' Plotting CDF - raw  <press any key>\n');
[f,xi] = ksdensity(ip_mh_rawChain_unified,'function','cdf');
plot(xi,f,'-b','linewidth',3)
title('Parameter Cumulative Distribution Function (raw chain)','fontname', 'Times', 'fontsize',20);
xlabel('Gravity (m/s^2)','fontname', 'Times', 'fontsize',20);
ylabel('CDF','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng sip_gravity_cdf_raw.png
pause;
clf;

%FILTERED
fprintf(1,' Plotting CDF - filtered  <press any key>\n');
[f,xi] = ksdensity(ip_mh_filtChain_unified,'function','cdf');
plot(xi,f,'-r','linewidth',3)
title('Parameter Cumulative Distribution Function (filtered chain)','fontname', 'Times', 'fontsize',20);
xlabel('Gravity (m/s^2)','fontname', 'Times', 'fontsize',20);
ylabel('CDF','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng sip_gravity_cdf_filt.png
pause;
clf;

% Autocorrelation plots ---------------------------------------------------
% RAW and FILTERED

%Let's use the autocorr function provided in
% http://sfb649.wiwi.hu-berlin.de/quantnet/index.php?p=show&id=900

fprintf(1,' Plotting autocorrelation  <press any key>\n');
nlags=10;
cd ../
[ACF_raw, lags, bounds] = autocorr_local(ip_mh_rawChain_unified, nlags, 0);
[ACF_filt, lags, bounds] = autocorr_local(ip_mh_filtChain_unified,nlags, 0);
plot(lags,ACF_raw,'bo-',lags,ACF_filt,'r*-','linewidth',3);
ylabel('Autocorrelation for \theta=g','fontname', 'Times', 'fontsize',20);
xlabel('Lag','fontname', 'Times', 'fontsize',20);
title('Parameter Autocorrelation','fontname', 'Times', 'fontsize',20);
grid on;
h=legend('Raw chain','Filtered chain','location','northeast');
set(h,'fontname', 'Times', 'fontsize',16);
set(gca,'FontSize',16);
print -dpng outputData/sip_gravity_autocorrelation_raw_filt.png
pause;
clf;

close;
%-------------------------------------------------------------------
% Autocorrelation plots ---------------------------------------------------
% If Matlab has Econometrics Toolbox, then you may use its own autocorr function
% 
% fprintf(1,' Plotting autocorrelation  <press any key>\n');
% nlags=10;
% [ACF_raw, lags, bounds] = autocorr(ip_mh_rawChain_unified, nlags, 0);
% [ACF_filt, lags, bounds] = autocorr(ip_mh_filtChain_unified,nlags, 0);
% plot(lags,ACF_raw,'bo-',lags,ACF_filt,'r*-','linewidth',3);
% ylabel('Autocorrelation for \theta=g','fontname', 'Times', 'fontsize',20);
% xlabel('Lag','fontname', 'Times', 'fontsize',20);
% title('Parameter Autocorrelation','fontname', 'Times', 'fontsize',20);
% grid on;
% h=legend('Raw chain','Filtered chain','location','northeast');
% set(h,'fontname', 'Times', 'fontsize',16);
% set(gca,'FontSize',16);
% print -dpng sip_gravity_autocorrelation_raw_filt.png
% pause;
% clf;
% cd matlab/
%-------------------------------------------------------------------



% Covariance matrix -------------------------------------------------------
fprintf(1, '\n Covariance matrix:  <press any key>\n');
cov_matrix_g = cov(ip_mh_rawChain_unified)



