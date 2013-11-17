clear all;
cd outputData/
ip_filt_chain
ip_raw_chain

% Chain plots -------------------------------------------------------------
% RAW
for i=1:size(ip_mh_filtChain_unified(),1)
 id_filt(i)=i;
end
 
%for i=1:size(ip_mh_rawChain_unified(),1)
% id_raw(i)=i;
% end
 
fprintf(1,' Plotting chains - filtered  <press any key>\n');
%plot(id_filt,ip_mh_filtChain_unified);
plot(id_filt,ip_mh_filtChain_unified(:,1),'b', id_filt,ip_mh_filtChain_unified(:,2),'g');
ylabel('Parameter values','fontname', 'Times', 'fontsize',20);
xlabel('Number of positions','fontname', 'Times', 'fontsize',20);
title('DRAM MCMC Chain Positions (filtered)','fontname', 'Times', 'fontsize',20);
h=legend('\theta_1','\theta_2','location','northeast');
axis([0 2000 -8 6]);
set(h,'fontname', 'Times', 'fontsize',16);
set(gca,'FontSize',16);
print -dpng simple_ip_chain_pos_filt.png
pause;
clf;

% Histogram plots ---------------------------------------------------------
% RAW
fprintf(1,' Plotting histogram - raw  <press any key>\n');
nbins=20;
cd ../
histyy(ip_mh_rawChain_unified(:,1),nbins,ip_mh_rawChain_unified(:,2),nbins, '\theta_1', '\theta_2')
h = findobj(gca,'Type','patch');
title('Parameters Histogram (raw chain, nbins=20)','fontname', 'Times', 'fontsize',20);
xlabel('Parameters','fontname', 'Times', 'fontsize',20);
grid on;
print -dpng outputData/simple_ip_hist_raw.png
pause;
clf;

%FILTERED 
fprintf(1,' Plotting histogram - filtered  <press any key>\n');
histyy(ip_mh_filtChain_unified(:,1),nbins,ip_mh_filtChain_unified(:,2),nbins, '\theta_1', '\theta_2')
h = findobj(gca,'Type','patch');
title('Parameters Histogram (filtered chain, nbins=20)','fontname', 'Times', 'fontsize',20);
xlabel('Parameters','fontname', 'Times', 'fontsize',20);
grid on;
print -dpng outputData/simple_ip_hist_filt.png
pause;
clf;
cd outputData/


% KDE plots ---------------------------------------------------------------
% RAW
fprintf(1,' Plotting KDE - raw  <press any key>\n');
[f,x] = ksdensity(ip_mh_rawChain_unified(:,1),'function','pdf');
[f2,x2] = ksdensity(ip_mh_rawChain_unified(:,2),'function','pdf');

x_p1=sort(ip_mh_rawChain_unified(:,1)); %analytical
f_p1=(exp(-(x_p1+1).*(x_p1+1)/8))/2/sqrt(2*pi);
x_p2=sort(ip_mh_rawChain_unified(:,1));
f_p2=(exp(-(x_p2-2).*(x_p2-2)/2))/sqrt(2*pi);

plot(x,f,'b',x2,f2,'g','linewidth',4);
hold;
plot(x_p1,f_p1,'--k',x_p2,f_p2,'-k','linewidth',2);

h=legend('\theta_1', '\theta_2', 'analytical (\theta_1)', 'analytical (\theta_2)', 'location', 'northwest');
title('Parameter Kernel Density Estimation (raw chain)','fontname', 'Times', 'fontsize',20);
xlabel('\theta_1 and \theta_2','fontname', 'Times', 'fontsize',20);
ylabel('KDE','fontname', 'Times', 'fontsize',20);
grid minor;
set(gca,'FontSize',16);
print -dpng simple_ip_kde_raw.png
pause;
clf;

%FILTERED
fprintf(1,' Plotting KDE - filtered  <press any key>\n');
[f,xi] = ksdensity(ip_mh_filtChain_unified(:,1),'function','pdf');
[f2,x2] = ksdensity(ip_mh_filtChain_unified(:,2),'function','pdf');

x_p1=sort(ip_mh_filtChain_unified(:,1)); %analytical
f_p1=(exp(-(x_p1+1).*(x_p1+1)/8))/2/sqrt(2*pi);
x_p2=sort(ip_mh_filtChain_unified(:,1));
f_p2=(exp(-(x_p2-2).*(x_p2-2)/2))/sqrt(2*pi);

plot(x,f,'b',x2,f2,'g','linewidth',4);
hold;
plot(x_p1,f_p1,'--k',x_p2,f_p2,'-k','linewidth',2);

h=legend('\theta_1', '\theta_2', 'analytical (\theta_1)', 'analytical (\theta_2)', 'location', 'northwest');
legend('\theta_1','\theta_2','location','northeast');
title('Parameter Kernel Density Estimation (filtered chain)','fontname', 'Times', 'fontsize',20);
xlabel('\theta_1 and \theta_2','fontname', 'Times', 'fontsize',20);
ylabel('KDE','fontname', 'Times', 'fontsize',20);
grid minor;
set(gca,'FontSize',16);
print -dpng simple_ip_kde_filt.png
pause;
clf;


% CDF plots ---------------------------------------------------------------
fprintf(1,' Plotting CDF - raw  <press any key>\n');
[f,xi] = ksdensity(ip_mh_rawChain_unified(:,1),'function','cdf');
[f2,x2] = ksdensity(ip_mh_rawChain_unified(:,2),'function','cdf');
plot(xi,f,'b',x2,f2,'g','linewidth',3);
h=legend('\theta_1','\theta_2','location','southeast');
title('Parameter Cumulative Distribution Function (raw chain)','fontname', 'Times', 'fontsize',20);
xlabel('\theta_1 and \theta_2','fontname', 'Times', 'fontsize',20);
ylabel('CDF','fontname', 'Times', 'fontsize',20);
grid minor;
set(gca,'FontSize',16);
print -dpng simple_ip_cdf_raw.png
pause;
clf;

%FILTERED
fprintf(1,' Plotting KCDF - filtered  <press any key>\n');
[f,xi] = ksdensity(ip_mh_filtChain_unified(:,1),'function','cdf');
[f2,x2] = ksdensity(ip_mh_filtChain_unified(:,2),'function','cdf');
plot(xi,f,'b',x2,f2,'g','linewidth',3);
h=legend('\theta_1','\theta_2','location','southeast');
title('Parameter Cumulative Distribution Function (filtered chain)','fontname', 'Times', 'fontsize',20);
xlabel('\theta_1 and \theta_2','fontname', 'Times', 'fontsize',20);
ylabel('CDF','fontname', 'Times', 'fontsize',20);
grid minor;
set(gca,'FontSize',16);
print -dpng simple_ip_cdf_filt.png
pause;
clf;

% Autocorrelation plots ---------------------------------------------------
% RAW and FILTERED

%Let's use the autocorr function provided in
% http://sfb649.wiwi.hu-berlin.de/quantnet/index.php?p=show&id=900

fprintf(1,' Plotting autocorrelation  <press any key>\n');
nlags=10;
[ACF_raw, lags, bounds] = autocorr(ip_mh_rawChain_unified(:,1), nlags, 0);
[ACF_filt, lags, bounds] = autocorr(ip_mh_filtChain_unified(:,1),nlags, 0);

[ACF_raw2, lags2, bounds2] = autocorr(ip_mh_rawChain_unified(:,2), nlags, 0);
[ACF_filt2, lags2, bounds2] = autocorr(ip_mh_filtChain_unified(:,2),nlags, 0);

plot(lags,ACF_raw,'b--*',lags,ACF_filt,'b*-',lags2,ACF_raw2,'g--*',lags2,ACF_filt2,'g*-','linewidth',3);

ylabel('Autocorrelation','fontname', 'Times', 'fontsize',20);
xlabel('Lag','fontname', 'Times', 'fontsize',20);
title('Parameter Autocorrelation','fontname', 'Times', 'fontsize',20);
h=legend('\theta_1, raw chain','\theta_1, filtered chain','\theta_2, raw chain','\theta_2, filtered chain','location','northeast');
set(h,'fontname', 'Times', 'fontsize',16);
axis([0 10 -.1 1]);
grid minor;
set(gca,'FontSize',16);
print -dpng simple_ip_autocorrelation_raw_filt.png
pause;
clf;

cd ..
close;


