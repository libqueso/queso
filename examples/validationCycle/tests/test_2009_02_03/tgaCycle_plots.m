cd outputData
tgaCalOutput_sub0
tgaCalExtraOutput_sub0
tgaValOutput_sub0
tgaValExtraOutput_sub0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
print -dpng unif_cal_prior_0.png
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(cycle_cal_ip_mc_filtChain_unifGkdePosits_sub0(1,:),cycle_cal_ip_mc_filtChain_unifGkdeValues_sub0(1,:),'b-','linewidth',3);
hold
plot(cycle_val_ip_mc_filtChain_unifGkdePosits_sub0(1,:),cycle_val_ip_mc_filtChain_unifGkdeValues_sub0(1,:),'r-','linewidth',3);
ylabel('Posterior marginal pdf','fontsize',20);
xlabel('A (min^{-1})','fontsize',20);
title('A: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','southwest');

plot(cycle_cal_ip_prior_0_grid,cycle_cal_ip_prior_0_values,'b*');

cal_left_vertical_line_x0 = ones(1,11)*minPrior0;
cal_left_vertical_line_y0 = [0 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(cal_left_vertical_line_x0,cal_left_vertical_line_y0,'b--','linewidth',1);

cal_right_vertical_line_x0 = ones(1,11)*maxPrior0;
cal_right_vertical_line_y0 = [0 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(cal_right_vertical_line_x0,cal_right_vertical_line_y0,'b--','linewidth',1);

print -dpng unif_cal_val_post_mpdf_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(cycle_cal_ip_mc_filtChain_unifGkdePosits_sub0(2,:),cycle_cal_ip_mc_filtChain_unifGkdeValues_sub0(2,:),'b-','linewidth',3);
hold
plot(cycle_val_ip_mc_filtChain_unifGkdePosits_sub0(2,:),cycle_val_ip_mc_filtChain_unifGkdeValues_sub0(2,:),'r-','linewidth',3);
ylabel('Pdf','fontsize',20);
xlabel('E (J/mol)','fontsize',20);
title('E: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');

plot(cycle_cal_ip_prior_1_grid,cycle_cal_ip_prior_1_values,'b*');

cal_left_vertical_line_x1 = ones(1,11)*minPrior1;
cal_left_vertical_line_y1 = [0 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(cal_left_vertical_line_x1,cal_left_vertical_line_y1,'b--','linewidth',1);

cal_right_vertical_line_x1 = ones(1,11)*maxPrior1;
cal_right_vertical_line_y1 = [0 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(cal_right_vertical_line_x1,cal_right_vertical_line_y1,'b--','linewidth',1);

print -dpng unif_cal_val_post_mpdf_1.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(cycle_cal_fp_mc_QoiSeq_unifGkdePosits_sub0(1,:),cycle_cal_fp_mc_QoiSeq_unifGkdeValues_sub0(1,:),'b-','linewidth',3);
hold
plot(cycle_val_fp_mc_QoiSeq_unifGkdePosits_sub0(1,:),cycle_val_fp_mc_QoiSeq_unifGkdeValues_sub0(1,:),'r-','linewidth',3);
ylabel('Pdf','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('QoI Pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng unif_cal_val_qoi_mpdf_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(cycle_cal_fp_unifQoiCdf_0_grid_sub0,100*cycle_cal_fp_unifQoiCdf_0_values_sub0,'b-','linewidth',3);
hold
plot(cycle_val_fp_unifQoiCdf_0_grid_sub0,100*cycle_val_fp_unifQoiCdf_0_values_sub0,'r-','linewidth',3);

a=axis;
axis([a(1) 1.0 a(3) a(4)]);
ylabel('Cdf (%)','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('QoI Cdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng unif_cal_val_qoi_Cdf_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(cycle_cal_fp_unifQoiCdf_0_grid_sub0,100*cycle_cal_fp_unifQoiCdf_0_values_sub0,'b-','linewidth',3);
hold
plot(cycle_val_fp_unifQoiCdf_0_grid_sub0,100*cycle_val_fp_unifQoiCdf_0_values_sub0,'r-','linewidth',3);

[m,n] = size(cycle_cal_fp_unifQoiCdf_0_grid_sub0);
epsilon = 0.06;
plot([0:0.01:1],ones(101,1)*100*(0.+epsilon/2),'k--','linewidth',1);
plot([0:0.01:1],ones(101,1)*100*(1.-epsilon/2),'k--','linewidth',1);

a=axis;
axis([a(1) 1.0 a(3) a(4)]);
ylabel('Cdf (%)','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('Full Picture','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng unif_cal_val_qoi_Cdf_0_model_confidence_a.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(cycle_cal_fp_unifQoiCdf_0_grid_sub0,100*cycle_cal_fp_unifQoiCdf_0_values_sub0,'b-','linewidth',3);
hold
plot(cycle_val_fp_unifQoiCdf_0_grid_sub0,100*cycle_val_fp_unifQoiCdf_0_values_sub0,'r-','linewidth',3);

[m,n] = size(cycle_cal_fp_unifQoiCdf_0_grid_sub0);
epsilon = 0.06;
plot([0:0.01:1],ones(101,1)*100*(0.+epsilon/2),'k--','linewidth',3);
plot([0:0.01:1],ones(101,1)*100*(1.-epsilon/2),'k--','linewidth',3);

a=axis;
axis([0.15 0.27 2 15]);
ylabel('Cdf (%)','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('Zoom','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng unif_cal_val_qoi_Cdf_0_model_confidence_b.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(cycle_cal_fp_unifQoiCdf_0_grid_sub0,100*cycle_cal_fp_unifQoiCdf_0_values_sub0,'b-','linewidth',3);
hold
plot(cycle_val_fp_unifQoiCdf_0_grid_sub0,100*cycle_val_fp_unifQoiCdf_0_values_sub0,'r-','linewidth',3);

[m,n] = size(cycle_cal_fp_unifQoiCdf_0_grid_sub0);
probabilityTresholdForFailure = 0.05;
plot([0:0.01:1],ones(101,1)*100*probabilityTresholdForFailure,'k--','linewidth',1);
massTresholdForFailure = 0.2;
plot(ones(101,1)*massTresholdForFailure,[0:1:100],'k--','linewidth',1);

a=axis;
axis([a(1) 1.0 a(3) a(4)]);
ylabel('Cdf (%)','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('Full Picture','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng unif_cal_val_qoi_Cdf_0_system_confidence_a.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(cycle_cal_fp_unifQoiCdf_0_grid_sub0,100*cycle_cal_fp_unifQoiCdf_0_values_sub0,'b-','linewidth',3);
hold
plot(cycle_val_fp_unifQoiCdf_0_grid_sub0,100*cycle_val_fp_unifQoiCdf_0_values_sub0,'r-','linewidth',3);

[m,n] = size(cycle_cal_fp_unifQoiCdf_0_grid_sub0);
probabilityTresholdForFailure = 0.05;
plot([0:0.01:1],ones(101,1)*100*probabilityTresholdForFailure,'k--','linewidth',3);
massTresholdForFailure = 0.2;
plot(ones(101,1)*massTresholdForFailure,[0:1:100],'k--','linewidth',3);

a=axis;
axis([0.15 0.25 2 8]);
ylabel('Cdf (%)','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('Zoom','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng unif_cal_val_qoi_Cdf_0_system_confidence_b.png
waitforbuttonpress;
clf

