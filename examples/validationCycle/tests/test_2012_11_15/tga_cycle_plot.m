cd outputData
%sfpOutput_sub0
%sfpExtraOutput_sub0
%cd ..
file_val_fp_qoi2;
file_cal_fp_qoi2;
file_cal_ip_raw;
file_val_ip_raw;

%Prior
%-------------------------------------------------------------------------

minPrior0 = 2.40e+11;
maxPrior0 = 2.80e+11;
numPrior0 = 200;
numHorizPts0 = 20;
deltaPrior0 = (maxPrior0 - minPrior0)/numPrior0;
cycle_cal_ip_prior_0_grid   = [minPrior0-numHorizPts0*deltaPrior0 : deltaPrior0 : maxPrior0+numHorizPts0*deltaPrior0];
cycle_cal_ip_prior_0_values = ones(1,numPrior0+2*numHorizPts0+1)./(maxPrior0-minPrior0);
cycle_cal_ip_prior_0_values(1,1:numHorizPts0)     = 0.;
cycle_cal_ip_prior_0_values(1,numPrior0+numHorizPts0+2:numPrior0+2*numHorizPts0+1) = 0.;
plot(cycle_cal_ip_prior_0_grid,cycle_cal_ip_prior_0_values,'b*');
hold

cal_left_vertical_line_x0 = ones(1,11)*minPrior0;
cal_left_vertical_line_y0 = [0 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(cal_left_vertical_line_x0,cal_left_vertical_line_y0,'b--','linewidth',1);

cal_right_vertical_line_x0 = ones(1,11)*maxPrior0;
cal_right_vertical_line_y0 = [0 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(cal_right_vertical_line_x0,cal_right_vertical_line_y0,'b--','linewidth',1);

ylabel('Prior marginal pdf','fontsize',20);
xlabel('A (min^{-1})','fontsize',20);
title('Parameter A: prior marginal pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng cal_prior_0_sub0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minPrior1 = 1.80e+5;
maxPrior1 = 2.20e+5;
numPrior1 = 200;
numHorizPts1 = 20;
deltaPrior1 = (maxPrior1 - minPrior1)/numPrior1;
cycle_cal_ip_prior_1_grid   = [minPrior1-numHorizPts1*deltaPrior1 : deltaPrior1 : maxPrior1+numHorizPts1*deltaPrior1];
cycle_cal_ip_prior_1_values = ones(1,numPrior1+2*numHorizPts1+1)./(maxPrior1-minPrior1);
cycle_cal_ip_prior_1_values(1,1:numHorizPts1)     = 0.;
cycle_cal_ip_prior_1_values(1,numPrior1+numHorizPts1+2:numPrior1+2*numHorizPts1+1) = 0.;
plot(cycle_cal_ip_prior_1_grid,cycle_cal_ip_prior_1_values,'b*');
hold

cal_left_vertical_line_x1 = ones(1,11)*minPrior1;
cal_left_vertical_line_y1 = [0 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(cal_left_vertical_line_x1,cal_left_vertical_line_y1,'b--','linewidth',1);

cal_right_vertical_line_x1 = ones(1,11)*maxPrior1;
cal_right_vertical_line_y1 = [0 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(cal_right_vertical_line_x1,cal_right_vertical_line_y1,'b--','linewidth',1);

ylabel('Prior marginal pdf','fontsize',20);
xlabel('E (J/mol)','fontsize',20);
title('Parameter E: prior marginal pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng unif_cal_prior_1.png
waitforbuttonpress;
clf




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
axis([0 1.0 0 5]);
set(gca,'fontsize',20);
ylabel('Pdf','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('QoI Pdf','fontsize',20);
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
a=axis;
axis([0 1.0 a(3) a(4)]);
grid minor;
set(gca,'fontsize',20);
title('QoI CDF');
ylabel('Cdf','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
legend('calibration',...
       'validation',...
       'location','northwest');
print -dpng fig3.png
waitforbuttonpress;
clf;

%CDF full picture

plot(xi,f,'r-',xiv,fv,'b-','linewidth',3);
hold
epsilon = 0.06;
plot([0:0.01:1],ones(101,1)*(0.+epsilon/2),'k--','linewidth',1);
plot([0:0.01:1],ones(101,1)*(1.-epsilon/2),'k--','linewidth',1);

a=axis;
axis([0 1.0 a(3) a(4)]);
ylabel('Cdf','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('Full Picture','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng cal_val_qoi_Cdf_0_model_confidence_a_sub0.png
waitforbuttonpress;
clf

%CDF Zoom

[f,xi] = ksdensity(cycle_cal_fp_mc_QoiSeq_unified,'function','cdf');
[fv,xiv] = ksdensity(cycle_val_fp_mc_QoiSeq_unified,'function','cdf');
plot(xi,f,'r-',xiv,fv,'b-','linewidth',3);
axis([0.15 0.27 .02 .15]);
ylabel('Cdf','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('Zoom','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng cal_val_qoi_Cdf_0_model_confidence_b_sub0.png
waitforbuttonpress;
clf


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