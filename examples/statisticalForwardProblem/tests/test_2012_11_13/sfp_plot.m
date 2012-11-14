cd outputData
%sfpOutput_sub0
%sfpExtraOutput_sub0
%cd ..
file_sfp_qoi

nlags=10;
[ACF, lags, bounds] = autocorr(fp_mc_QoiSeq_unified, nlags, 0);
plot(lags,ACF,'mo-','linewidth',3);
title('QoI autocorrelation');
print -dpng fig1.png
waitforbuttonpress;
clf;

[f,xi] = ksdensity(fp_mc_QoiSeq_unified,'function','pdf');
plot(xi,f,'-m','linewidth',3)
title('QoI KDE');
print -dpng fig2.png
waitforbuttonpress;
clf;

[f,xi] = ksdensity(fp_mc_QoiSeq_unified,'function','cdf');
plot(xi,f,'-m','linewidth',3)
title('QoI cdf');
print -dpng fig3.png
waitforbuttonpress;
clf;

cd ..