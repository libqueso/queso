clear all;
cd outputData/
fp_p_seq
fp_q_seq

% Chain plots -------------------------------------------------------------
% Parameter
for i=1:size(fp_mc_ParamSeq_unified(),1)
 id_p(i)=i;
end

fprintf(1,' Plotting chains - parameter  <press any key>\n');

plot(id_p,fp_mc_ParamSeq_unified(:,1),'b', id_p,fp_mc_ParamSeq_unified(:,2),'g');
ylabel('Parameter values','fontname', 'Times', 'fontsize',20);
xlabel('Number of positions','fontname', 'Times', 'fontsize',20);
title('DRAM MCMC Chain Positions (filtered)','fontname', 'Times', 'fontsize',20);
h=legend('\theta_1','\theta_2','location','northeast');
axis([0 2000 -8 6]);
set(h,'fontname', 'Times', 'fontsize',16);
set(gca,'FontSize',16);
print -dpng simple_fp_chain_pos_param.png
pause;
clf;

%QOI
fprintf(1,' Plotting chains - QoI  <press any key>\n');

plot(fp_mc_QoiSeq_unified,'m');
ylabel('QoI values','fontname', 'Times', 'fontsize',20);
xlabel('Number of positions','fontname', 'Times', 'fontsize',20);
title('DRAM MCMC Chain Positions - QoI','fontname', 'Times', 'fontsize',20);
set(gca,'FontSize',16);
print -dpng simple_fp_chain_pos_qoi.png
pause;
clf;

% Histogram plots ---------------------------------------------------------
fprintf(1,' Plotting histogram - raw  <press any key>\n');
nbins=20;
hist(fp_mc_QoiSeq_unified,nbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','m','EdgeColor','k');%
grid on;
xlabel('QoI=\theta_1+\theta_2','fontname', 'Times', 'fontsize',20);
ylabel('Frequency','fontname', 'Times', 'fontsize',20);
print -dpng simple_fp_hist_qoi.png
pause;
clf;

% KDE plots ---------------------------------------------------------------
[fi,xi] = ksdensity(fp_mc_QoiSeq_unified,'function','pdf');
x=sort(fp_mc_QoiSeq_unified);
mu=1;
sigma2=5;
f=(exp(-(x-mu).*(x-mu)/sigma2/2))/sqrt(2*pi*sigma2);
plot(xi,fi,'-m','linewidth',4);
hold;
plot(x,f,'--k','linewidth',2);
h=legend('QoI = \theta_1+\theta_2', 'analytical', 'location', 'northwest');
title('QoI Kernel Density Estimation ','fontname', 'Times', 'fontsize',20);
xlabel('QoI','fontname', 'Times', 'fontsize',20);
ylabel('KDE','fontname', 'Times', 'fontsize',20);
grid minor;
set(gca,'FontSize',16);
print -dpng simple_fp_kde_qoi.png
hold off;
pause;
clf;


% CDF plots ---------------------------------------------------------------
[f,xi] = ksdensity(fp_mc_QoiSeq_unified,'function','cdf');
plot(xi,f,'-m','linewidth',3)
title('QoI Cumulative Distribution Function','fontname', 'Times', 'fontsize',20);
xlabel('QoI=\theta_1 + \theta_2','fontname', 'Times', 'fontsize',20);
ylabel('CDF','fontname', 'Times', 'fontsize',20);
grid minor;
set(gca,'FontSize',16);
print -dpng simple_fp_cdf_qoi.png
pause;
clf;

% Autocorrelation plots ---------------------------------------------------
% RAW and FILTERED

%Let's use the autocorr function provided in
% http://sfb649.wiwi.hu-berlin.de/quantnet/index.php?p=show&id=900

fprintf(1,' Plotting autocorrelation  <press any key>\n');
nlags=10;
[ACF, lags, bounds] = autocorr(fp_mc_QoiSeq_unified, nlags, 0);

plot(lags,ACF,'m--*', 'linewidth',3);

ylabel('Autocorrelation','fontname', 'Times', 'fontsize',20);
xlabel('Lag','fontname', 'Times', 'fontsize',20);
title('QoI Autocorrelation','fontname', 'Times', 'fontsize',20);
axis([0 10 -.1 1]);
grid minor;
set(gca,'FontSize',16);
print -dpng simple_fp_autocorrelation_qoi.png
pause;
clf;

cd ..
close;


