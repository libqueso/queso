cd outputData
%sfpOutput_sub0
%sfpExtraOutput_sub0
%cd ..
file_val_fp_qoi2;
file_cal_fp_qoi2;
file_cal_ip_raw;
file_val_ip_raw;
%FORWARD
%-------------------------------------------------------------------------
% Autocorrelation
nlags=10;
[ACF, lags, bounds] = autocorr(cycle_cal_fp_mc_QoiSeq_unified, nlags, 0);
[ACFv, lagsv, boundsv] = autocorr(cycle_val_fp_mc_QoiSeq_unified, nlags, 0);
plot(lags,ACF,'ro-',lagsv,ACFv,'b*-','linewidth',3);
grid minor;
set(gca,'fontsize',20);
title('QoI autocorrelation');
legend('calibration',...
       'validation',...
       'location','northwest');
print -dpng fig1.png
waitforbuttonpress;
clf;

%KDE
[f,xi] = ksdensity(cycle_cal_fp_mc_QoiSeq_unified,'function','pdf');
[fv,xiv] = ksdensity(cycle_val_fp_mc_QoiSeq_unified,'function','pdf');
plot(xi,f,'r-',xiv,fv,'b-','linewidth',3);
grid minor;
set(gca,'fontsize',20);
title('QoI KDE');
legend('calibration',...
       'validation',...
       'location','northwest');
print -dpng fig2.png
waitforbuttonpress;
clf;

clear f;
clear xi;
clear fv;
clear xiv;
%CDF
[f,xi] = ksdensity(cycle_cal_fp_mc_QoiSeq_unified,'function','cdf');
[fv,xiv] = ksdensity(cycle_val_fp_mc_QoiSeq_unified,'function','cdf');
plot(xi,f,'r-',xiv,fv,'b-','linewidth',3);
grid minor;
set(gca,'fontsize',20);
title('QoI CDF');
legend('calibration',...
       'validation',...
       'location','northwest');
print -dpng fig3.png
waitforbuttonpress;
clf;

%INVERSE
%------------------------------------------------------------------------
% Parameter 1 ------------
% KDE plots 
[f1,x1] = ksdensity(cycle_cal_ip_mh_rawChain_unified(1,:),'function','pdf');
[f2,x2] = ksdensity(cycle_val_ip_mh_rawChain_unified(1,:),'function','pdf');
plot(x1,f1,'-b',x2,f2,'-r','linewidth',3)
title('Posterior marginal pdfs of parameter 1','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('calibration',...
       'validation',...
       'location','northwest');
print -dpng fig4.png
waitforbuttonpress;
clf;


% CDF plots 
[f1,x1] = ksdensity(cycle_cal_ip_mh_rawChain_unified(1,:),'function','cdf');
[f2,x2] = ksdensity(cycle_val_ip_mh_rawChain_unified(1,:),'function','cdf');
plot(x1,f1,'-b',x2,f2,'-r','linewidth',3)
title('Cdfs of parameter 1','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('calibration',...
       'validation',...
       'location','northwest');
print -dpng fig5.png
waitforbuttonpress;
clf;
% Parameter 2 ------------
% KDE plots 
[f1,x1] = ksdensity(cycle_cal_ip_mh_rawChain_unified(2,:),'function','pdf');
[f2,x2] = ksdensity(cycle_val_ip_mh_rawChain_unified(2,:),'function','pdf');
plot(x1,f1,'-b',x2,f2,'-r','linewidth',3)
title('Posterior marginal pdfs of parameter 2','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('calibration',...
       'validation',...
       'location','northwest');
print -dpng fig6.png
waitforbuttonpress;
clf;


% CDF plots 
[f1,x1] = ksdensity(cycle_cal_ip_mh_rawChain_unified(2,:),'function','cdf');
[f2,x2] = ksdensity(cycle_val_ip_mh_rawChain_unified(2,:),'function','cdf');
plot(x1,f1,'-b',x2,f2,'-r','linewidth',3)
title('Cdfs of parameter 2','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('calibration',...
       'validation',...
       'location','northwest');
print -dpng fig7.png
waitforbuttonpress;
clf;

cd ..