cd outputData
tgaCalOutput_sub0
%tgaCalExtraOutput_sub0
tgaValOutput_sub0
%tgaValExtraOutput_sub0
tgaCalOutput_sub1
%tgaCalExtraOutput_sub1
tgaValOutput_sub1
%tgaValExtraOutput_sub1

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

cal_left_vertical_line_x0 = ones(1,11)*minPrior0;
cal_left_vertical_line_y0 = [0 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)];

cal_right_vertical_line_x0 = ones(1,11)*maxPrior0;
cal_right_vertical_line_y0 = [0 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)];

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

cal_left_vertical_line_x1 = ones(1,11)*minPrior1;
cal_left_vertical_line_y1 = [0 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)];

cal_right_vertical_line_x1 = ones(1,11)*maxPrior1;
cal_right_vertical_line_y1 = [0 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(cycle_cal_ip_mc_filtChain_unifGkdePosits_sub0(1,:),cycle_cal_ip_mc_filtChain_unifGkdeValues_sub0(1,:),'bo','linewidth',2);
hold
plot(cycle_val_ip_mc_filtChain_unifGkdePosits_sub0(1,:),cycle_val_ip_mc_filtChain_unifGkdeValues_sub0(1,:),'ro','linewidth',2);
plot(cycle_cal_ip_mc_filtChain_GkdePosits_sub0(1,:),cycle_cal_ip_mc_filtChain_GkdeValues_sub0(1,:),'b-','linewidth',2);
plot(cycle_val_ip_mc_filtChain_GkdePosits_sub0(1,:),cycle_val_ip_mc_filtChain_GkdeValues_sub0(1,:),'r-','linewidth',2);
plot(cycle_cal_ip_mc_filtChain_GkdePosits_sub1(1,:),cycle_cal_ip_mc_filtChain_GkdeValues_sub1(1,:),'b--','linewidth',2);
plot(cycle_val_ip_mc_filtChain_GkdePosits_sub1(1,:),cycle_val_ip_mc_filtChain_GkdeValues_sub1(1,:),'r--','linewidth',2);
ylabel('Posterior marginal pdf','fontsize',20);
xlabel('A (min^{-1})','fontsize',20);
title('A: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Unif Calibration',...
       'Unif Validation',...
       'Sub0 Calibration',...
       'Sub0 Validation',...
       'Sub1 Calibration',...
       'Sub1 Validation',...
       'location','southwest');

plot(cycle_cal_ip_prior_0_grid,cycle_cal_ip_prior_0_values,'b*');

cal_left_vertical_line_x0 = ones(1,11)*minPrior0;
cal_left_vertical_line_y0 = [0 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(cal_left_vertical_line_x0,cal_left_vertical_line_y0,'b--','linewidth',1);

cal_right_vertical_line_x0 = ones(1,11)*maxPrior0;
cal_right_vertical_line_y0 = [0 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : cycle_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(cal_right_vertical_line_x0,cal_right_vertical_line_y0,'b--','linewidth',1);

print -dpng compare_cal_val_post_mpdf_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(cycle_cal_ip_mc_filtChain_unifGkdePosits_sub0(2,:),cycle_cal_ip_mc_filtChain_unifGkdeValues_sub0(2,:),'bo','linewidth',2);
hold
plot(cycle_val_ip_mc_filtChain_unifGkdePosits_sub0(2,:),cycle_val_ip_mc_filtChain_unifGkdeValues_sub0(2,:),'ro','linewidth',2);
plot(cycle_cal_ip_mc_filtChain_GkdePosits_sub0(2,:),cycle_cal_ip_mc_filtChain_GkdeValues_sub0(2,:),'b-','linewidth',2);
plot(cycle_val_ip_mc_filtChain_GkdePosits_sub0(2,:),cycle_val_ip_mc_filtChain_GkdeValues_sub0(2,:),'r-','linewidth',2);
plot(cycle_cal_ip_mc_filtChain_GkdePosits_sub1(2,:),cycle_cal_ip_mc_filtChain_GkdeValues_sub1(2,:),'b--','linewidth',2);
plot(cycle_val_ip_mc_filtChain_GkdePosits_sub1(2,:),cycle_val_ip_mc_filtChain_GkdeValues_sub1(2,:),'r--','linewidth',2);
ylabel('Pdf','fontsize',20);
xlabel('E (J/mol)','fontsize',20);
title('E: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Unif Calibration',...
       'Unif Validation',...
       'Sub0 Calibration',...
       'Sub0 Validation',...
       'Sub1 Calibration',...
       'Sub1 Validation',...
       'location','northwest');

plot(cycle_cal_ip_prior_1_grid,cycle_cal_ip_prior_1_values,'b*');

cal_left_vertical_line_x1 = ones(1,11)*minPrior1;
cal_left_vertical_line_y1 = [0 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(cal_left_vertical_line_x1,cal_left_vertical_line_y1,'b--','linewidth',1);

cal_right_vertical_line_x1 = ones(1,11)*maxPrior1;
cal_right_vertical_line_y1 = [0 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : cycle_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(cal_right_vertical_line_x1,cal_right_vertical_line_y1,'b--','linewidth',1);

print -dpng compare_cal_val_post_mpdf_1.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(cycle_cal_fp_mc_QoiSeq_unifGkdePosits_sub0(1,:),cycle_cal_fp_mc_QoiSeq_unifGkdeValues_sub0(1,:),'bo','linewidth',2);
hold
plot(cycle_val_fp_mc_QoiSeq_unifGkdePosits_sub0(1,:),cycle_val_fp_mc_QoiSeq_unifGkdeValues_sub0(1,:),'ro','linewidth',2);
plot(cycle_cal_fp_mc_QoiSeq_GkdePosits_sub0(1,:),cycle_cal_fp_mc_QoiSeq_GkdeValues_sub0(1,:),'b-','linewidth',2);
plot(cycle_val_fp_mc_QoiSeq_GkdePosits_sub0(1,:),cycle_val_fp_mc_QoiSeq_GkdeValues_sub0(1,:),'r-','linewidth',2);
plot(cycle_cal_fp_mc_QoiSeq_GkdePosits_sub1(1,:),cycle_cal_fp_mc_QoiSeq_GkdeValues_sub1(1,:),'b--','linewidth',2);
plot(cycle_val_fp_mc_QoiSeq_GkdePosits_sub1(1,:),cycle_val_fp_mc_QoiSeq_GkdeValues_sub1(1,:),'r--','linewidth',2);
ylabel('Pdf','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('QoI Pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Unif Calibration',...
       'Unif Validation',...
       'Sub0 Calibration',...
       'Sub0 Validation',...
       'Sub1 Calibration',...
       'Sub1 Validation',...
       'location','northwest');
print -dpng compare_cal_val_qoi_mpdf_0.png
waitforbuttonpress;
clf
