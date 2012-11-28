cd outputData
file_sfp_qoi

% Plotting QoI autocorrelation -------------------------------------------
nlags=10;
[ACF, lags, bounds] = autocorr(fp_mc_QoiSeq_unified, nlags, 0);
plot(lags,ACF,'mo-','linewidth',3);
grid minor;
set(gca,'fontsize',20);
title('QoI autocorrelation');
print -dpng QoI_autocorrelation.png
waitforbuttonpress;
clf;

% Plotting QoI marginal posterior pdf ------------------------------------
[f,xi] = ksdensity(fp_mc_QoiSeq_unified,'function','pdf');
plot(xi,f,'-m','linewidth',3);
set(gca,'fontsize',20);
grid minor;
title('QoI KDE');
print -dpng QoI_PDF.png
waitforbuttonpress;
clf;

% Plotting QoI marginal posterior cdf ------------------------------------
[f,xi] = ksdensity(fp_mc_QoiSeq_unified,'function','cdf');
plot(xi,f,'-m','linewidth',3);
set(gca,'fontsize',20);
grid minor;
title('QoI cdf');
print -dpng QoI_CDF.png
waitforbuttonpress;
clf;

cd ..