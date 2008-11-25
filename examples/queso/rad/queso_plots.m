s1Output
s1ExtraOutput
s2Output
s2ExtraOutput

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minPrior0 = 1.; %2.40e+11;
maxPrior0 = 1.2; %2.80e+11;
numPrior0 = 200;
numHorizPts0 = 20;
deltaPrior0 = (maxPrior0 - minPrior0)/numPrior0;
rad_cal_ip_prior_0_grid   = [minPrior0-numHorizPts0*deltaPrior0 : deltaPrior0 : maxPrior0+numHorizPts0*deltaPrior0];
rad_cal_ip_prior_0_values = ones(1,numPrior0+2*numHorizPts0+1)./(maxPrior0-minPrior0);
rad_cal_ip_prior_0_values(1,1:numHorizPts0)     = 0.;
rad_cal_ip_prior_0_values(1,numPrior0+numHorizPts0+2:numPrior0+2*numHorizPts0+1) = 0.;
plot(rad_cal_ip_prior_0_grid,rad_cal_ip_prior_0_values,'b*');
hold

rad_cal_left_vertical_line_x0 = ones(1,11)*minPrior0;
rad_cal_left_vertical_line_y0 = [0 : rad_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : rad_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(rad_cal_left_vertical_line_x0,rad_cal_left_vertical_line_y0,'b--','linewidth',1);

rad_cal_right_vertical_line_x0 = ones(1,11)*maxPrior0;
rad_cal_right_vertical_line_y0 = [0 : rad_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : rad_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(rad_cal_right_vertical_line_x0,rad_cal_right_vertical_line_y0,'b--','linewidth',1);

ylabel('Prior marginal pdf','fontsize',20);
xlabel('CP1 [(atm m)^{-1}]','fontsize',20);
title('Parameter CP1: prior marginal pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng talk_rad_cal_prior_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minPrior1 = .004; %1.80e+5;
maxPrior1 = .006; %2.20e+5;
numPrior1 = 200;
numHorizPts1 = 20;
deltaPrior1 = (maxPrior1 - minPrior1)/numPrior1;
rad_cal_ip_prior_1_grid   = [minPrior1-numHorizPts1*deltaPrior1 : deltaPrior1 : maxPrior1+numHorizPts1*deltaPrior1];
rad_cal_ip_prior_1_values = ones(1,numPrior1+2*numHorizPts1+1)./(maxPrior1-minPrior1);
rad_cal_ip_prior_1_values(1,1:numHorizPts1)     = 0.;
rad_cal_ip_prior_1_values(1,numPrior1+numHorizPts1+2:numPrior1+2*numHorizPts1+1) = 0.;
plot(rad_cal_ip_prior_1_grid,rad_cal_ip_prior_1_values,'b*');
hold

rad_cal_left_vertical_line_x1 = ones(1,11)*minPrior1;
rad_cal_left_vertical_line_y1 = [0 : rad_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : rad_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(rad_cal_left_vertical_line_x1,rad_cal_left_vertical_line_y1,'b--','linewidth',1);

rad_cal_right_vertical_line_x1 = ones(1,11)*maxPrior1;
rad_cal_right_vertical_line_y1 = [0 : rad_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : rad_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(rad_cal_right_vertical_line_x1,rad_cal_right_vertical_line_y1,'b--','linewidth',1);

ylabel('Prior marginal pdf','fontsize',20);
xlabel('CP2 [m^{-1}]','fontsize',20);
title('Parameter CP2: prior marginal pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng talk_rad_cal_prior_1.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minPrior2 = 1.6e-10; %1.80e+5;
maxPrior2 = 2e-10; %2.20e+5;
numPrior2 = 200;
numHorizPts2 = 20;
deltaPrior2 = (maxPrior2 - minPrior2)/numPrior2;
rad_cal_ip_prior_2_grid   = [minPrior2-numHorizPts2*deltaPrior2 : deltaPrior2 : maxPrior2+numHorizPts2*deltaPrior2];
rad_cal_ip_prior_2_values = ones(1,numPrior2+2*numHorizPts2+1)./(maxPrior2-minPrior2);
rad_cal_ip_prior_2_values(1,1:numHorizPts2)     = 0.;
rad_cal_ip_prior_2_values(1,numPrior2+numHorizPts2+2:numPrior2+2*numHorizPts2+1) = 0.;
plot(rad_cal_ip_prior_2_grid,rad_cal_ip_prior_2_values,'b*');
hold

rad_cal_left_vertical_line_x2 = ones(1,11)*minPrior2;
rad_cal_left_vertical_line_y2 = [0 : rad_cal_ip_prior_2_values(1,numHorizPts2+1)/10 : rad_cal_ip_prior_2_values(1,numHorizPts2+1)];
plot(rad_cal_left_vertical_line_x2,rad_cal_left_vertical_line_y2,'b--','linewidth',1);

rad_cal_right_vertical_line_x2 = ones(1,11)*maxPrior2;
rad_cal_right_vertical_line_y2 = [0 : rad_cal_ip_prior_2_values(1,numHorizPts2+1)/10 : rad_cal_ip_prior_2_values(1,numHorizPts2+1)];
plot(rad_cal_right_vertical_line_x2,rad_cal_right_vertical_line_y2,'b--','linewidth',1);

ylabel('Prior marginal pdf','fontsize',20);
xlabel('CT1 [K^{-2}]','fontsize',20);
title('Parameter CT1: prior marginal pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng talk_rad_cal_prior_2.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minPrior3 = -5.e-7; %1.80e+5;
maxPrior3 = -4.e-7; %2.20e+5;
numPrior3 = 200;
numHorizPts3 = 20;
deltaPrior3 = (maxPrior3 - minPrior3)/numPrior3;
rad_cal_ip_prior_3_grid   = [minPrior3-numHorizPts3*deltaPrior3 : deltaPrior3 : maxPrior3+numHorizPts3*deltaPrior3];
rad_cal_ip_prior_3_values = ones(1,numPrior3+2*numHorizPts3+1)./(maxPrior3-minPrior3);
rad_cal_ip_prior_3_values(1,1:numHorizPts3)     = 0.;
rad_cal_ip_prior_3_values(1,numPrior3+numHorizPts3+2:numPrior3+2*numHorizPts3+1) = 0.;
plot(rad_cal_ip_prior_3_grid,rad_cal_ip_prior_3_values,'b*');
hold

rad_cal_left_vertical_line_x3 = ones(1,11)*minPrior3;
rad_cal_left_vertical_line_y3 = [0 : rad_cal_ip_prior_3_values(1,numHorizPts3+1)/10 : rad_cal_ip_prior_3_values(1,numHorizPts3+1)];
plot(rad_cal_left_vertical_line_x3,rad_cal_left_vertical_line_y3,'b--','linewidth',1);

rad_cal_right_vertical_line_x3 = ones(1,11)*maxPrior3;
rad_cal_right_vertical_line_y3 = [0 : rad_cal_ip_prior_3_values(1,numHorizPts3+1)/10 : rad_cal_ip_prior_3_values(1,numHorizPts3+1)];
plot(rad_cal_right_vertical_line_x3,rad_cal_right_vertical_line_y3,'b--','linewidth',1);

ylabel('Prior marginal pdf','fontsize',20);
xlabel('CT2 ([K^{-1}]','fontsize',20);
title('Parameter CT2: prior marginal pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng talk_rad_cal_prior_3.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minPrior4 = .004; %1.80e+5;
maxPrior4 = .006; %2.20e+5;
numPrior4 = 200;
numHorizPts4 = 20;
deltaPrior4 = (maxPrior4 - minPrior4)/numPrior4;
rad_cal_ip_prior_4_grid   = [minPrior4-numHorizPts4*deltaPrior4 : deltaPrior4 : maxPrior4+numHorizPts4*deltaPrior4];
rad_cal_ip_prior_4_values = ones(1,numPrior4+2*numHorizPts4+1)./(maxPrior4-minPrior4);
rad_cal_ip_prior_4_values(1,1:numHorizPts4)     = 0.;
rad_cal_ip_prior_4_values(1,numPrior4+numHorizPts4+2:numPrior4+2*numHorizPts4+1) = 0.;
plot(rad_cal_ip_prior_4_grid,rad_cal_ip_prior_4_values,'b*');
hold

rad_cal_left_vertical_line_x4 = ones(1,11)*minPrior4;
rad_cal_left_vertical_line_y4 = [0 : rad_cal_ip_prior_4_values(1,numHorizPts4+1)/10 : rad_cal_ip_prior_4_values(1,numHorizPts4+1)];
plot(rad_cal_left_vertical_line_x4,rad_cal_left_vertical_line_y4,'b--','linewidth',1);

rad_cal_right_vertical_line_x4 = ones(1,11)*maxPrior4;
rad_cal_right_vertical_line_y4 = [0 : rad_cal_ip_prior_4_values(1,numHorizPts4+1)/10 : rad_cal_ip_prior_4_values(1,numHorizPts4+1)];
plot(rad_cal_right_vertical_line_x4,rad_cal_right_vertical_line_y4,'b--','linewidth',1);

ylabel('Prior marginal pdf','fontsize',20);
xlabel('CT3 [non-dimensional]','fontsize',20);
title('Parameter CT3: prior marginal pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng talk_rad_cal_prior_4.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(rad_cal_ip_mdf_0_grid,rad_cal_ip_mdf_0_values,'b+');
plot(rad_cal_ip_mc_filteredChain_kdeEvalPositions(1,:),rad_cal_ip_mc_filteredChain_gaussianKdeDensities(1,:),'b-','linewidth',3);
hold
%%plot(rad_val_ip_mdf_0_grid,rad_val_ip_mdf_0_values,'r+');
plot(rad_val_ip_mc_filteredChain_kdeEvalPositions(1,:),rad_val_ip_mc_filteredChain_gaussianKdeDensities(1,:),'r-','linewidth',3);
ylabel('Posterior marginal pdf','fontsize',20);
xlabel('CP1 [(atm m){-1}]','fontsize',20);
title('CP1: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','southwest');

plot(rad_cal_ip_prior_0_grid,rad_cal_ip_prior_0_values,'b*');

rad_cal_left_vertical_line_x0 = ones(1,11)*minPrior0;
rad_cal_left_vertical_line_y0 = [0 : rad_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : rad_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(rad_cal_left_vertical_line_x0,rad_cal_left_vertical_line_y0,'b--','linewidth',1);

rad_cal_right_vertical_line_x0 = ones(1,11)*maxPrior0;
rad_cal_right_vertical_line_y0 = [0 : rad_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : rad_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(rad_cal_right_vertical_line_x0,rad_cal_right_vertical_line_y0,'b--','linewidth',1);

print -dpng talk_rad_cal_val_post_mpdf_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(rad_cal_ip_mdf_1_grid,rad_cal_ip_mdf_1_values,'b+');
plot(rad_cal_ip_mc_filteredChain_kdeEvalPositions(2,:),rad_cal_ip_mc_filteredChain_gaussianKdeDensities(2,:),'b-','linewidth',3);
hold
%%plot(rad_val_ip_mdf_1_grid,rad_val_ip_mdf_1_values,'r+');
plot(rad_val_ip_mc_filteredChain_kdeEvalPositions(2,:),rad_val_ip_mc_filteredChain_gaussianKdeDensities(2,:),'r-','linewidth',3);
ylabel('Pdf','fontsize',20);
xlabel('CP2 [m^{-1}]','fontsize',20);
title('CP2: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');

plot(rad_cal_ip_prior_1_grid,rad_cal_ip_prior_1_values,'b*');

rad_cal_left_vertical_line_x1 = ones(1,11)*minPrior1;
rad_cal_left_vertical_line_y1 = [0 : rad_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : rad_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(rad_cal_left_vertical_line_x1,rad_cal_left_vertical_line_y1,'b--','linewidth',1);

rad_cal_right_vertical_line_x1 = ones(1,11)*maxPrior1;
rad_cal_right_vertical_line_y1 = [0 : rad_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : rad_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(rad_cal_right_vertical_line_x1,rad_cal_right_vertical_line_y1,'b--','linewidth',1);

print -dpng talk_rad_cal_val_post_mpdf_1.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(rad_cal_ip_mdf_0_grid,rad_cal_ip_mdf_0_values,'b+');
plot(rad_cal_ip_mc_filteredChain_kdeEvalPositions(3,:),rad_cal_ip_mc_filteredChain_gaussianKdeDensities(3,:),'b-','linewidth',3);
hold
%%plot(rad_val_ip_mdf_0_grid,rad_val_ip_mdf_0_values,'r+');
plot(rad_val_ip_mc_filteredChain_kdeEvalPositions(3,:),rad_val_ip_mc_filteredChain_gaussianKdeDensities(3,:),'r-','linewidth',3);
ylabel('Posterior marginal pdf','fontsize',20);
xlabel('CT1 [K^{-2}]','fontsize',20);
title('CT1: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','southwest');

plot(rad_cal_ip_prior_2_grid,rad_cal_ip_prior_2_values,'b*');

rad_cal_left_vertical_line_x2 = ones(1,11)*minPrior2;
rad_cal_left_vertical_line_y2 = [0 : rad_cal_ip_prior_2_values(1,numHorizPts2+1)/10 : rad_cal_ip_prior_2_values(1,numHorizPts2+1)];
plot(rad_cal_left_vertical_line_x2,rad_cal_left_vertical_line_y2,'b--','linewidth',1);

rad_cal_right_vertical_line_x2 = ones(1,11)*maxPrior2;
rad_cal_right_vertical_line_y2 = [0 : rad_cal_ip_prior_2_values(1,numHorizPts2+1)/10 : rad_cal_ip_prior_2_values(1,numHorizPts2+1)];
plot(rad_cal_right_vertical_line_x2,rad_cal_right_vertical_line_y2,'b--','linewidth',1);

print -dpng talk_rad_cal_val_post_mpdf_2.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(rad_cal_ip_mdf_0_grid,rad_cal_ip_mdf_0_values,'b+');
plot(rad_cal_ip_mc_filteredChain_kdeEvalPositions(4,:),rad_cal_ip_mc_filteredChain_gaussianKdeDensities(4,:),'b-','linewidth',3);
hold
%%plot(rad_val_ip_mdf_0_grid,rad_val_ip_mdf_0_values,'r+');
plot(rad_val_ip_mc_filteredChain_kdeEvalPositions(4,:),rad_val_ip_mc_filteredChain_gaussianKdeDensities(4,:),'r-','linewidth',3);
ylabel('Posterior marginal pdf','fontsize',20);
xlabel('CT2 [K^{-1}]','fontsize',20);
title('CT2: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','southwest');

plot(rad_cal_ip_prior_3_grid,rad_cal_ip_prior_3_values,'b*');

rad_cal_left_vertical_line_x3 = ones(1,11)*minPrior3;
rad_cal_left_vertical_line_y3 = [0 : rad_cal_ip_prior_3_values(1,numHorizPts3+1)/10 : rad_cal_ip_prior_3_values(1,numHorizPts3+1)];
plot(rad_cal_left_vertical_line_x3,rad_cal_left_vertical_line_y3,'b--','linewidth',1);

rad_cal_right_vertical_line_x3 = ones(1,11)*maxPrior3;
rad_cal_right_vertical_line_y3 = [0 : rad_cal_ip_prior_3_values(1,numHorizPts3+1)/10 : rad_cal_ip_prior_3_values(1,numHorizPts3+1)];
plot(rad_cal_right_vertical_line_x3,rad_cal_right_vertical_line_y3,'b--','linewidth',1);

print -dpng talk_rad_cal_val_post_mpdf_3.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(rad_cal_ip_mdf_0_grid,rad_cal_ip_mdf_0_values,'b+');
plot(rad_cal_ip_mc_filteredChain_kdeEvalPositions(5,:),rad_cal_ip_mc_filteredChain_gaussianKdeDensities(5,:),'b-','linewidth',3);
hold
%%plot(rad_val_ip_mdf_0_grid,rad_val_ip_mdf_0_values,'r+');
plot(rad_val_ip_mc_filteredChain_kdeEvalPositions(5,:),rad_val_ip_mc_filteredChain_gaussianKdeDensities(5,:),'r-','linewidth',3);
ylabel('Posterior marginal pdf','fontsize',20);
xlabel('CT3 [non-dimensional]','fontsize',20);
title('CT3: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','southwest');

plot(rad_cal_ip_prior_4_grid,rad_cal_ip_prior_4_values,'b*');

rad_cal_left_vertical_line_x4 = ones(1,11)*minPrior4;
rad_cal_left_vertical_line_y4 = [0 : rad_cal_ip_prior_4_values(1,numHorizPts4+1)/10 : rad_cal_ip_prior_4_values(1,numHorizPts4+1)];
plot(rad_cal_left_vertical_line_x4,rad_cal_left_vertical_line_y4,'b--','linewidth',1);

rad_cal_right_vertical_line_x4 = ones(1,11)*maxPrior4;
rad_cal_right_vertical_line_y4 = [0 : rad_cal_ip_prior_4_values(1,numHorizPts4+1)/10 : rad_cal_ip_prior_4_values(1,numHorizPts4+1)];
plot(rad_cal_right_vertical_line_x4,rad_cal_right_vertical_line_y4,'b--','linewidth',1);

print -dpng talk_rad_cal_val_post_mpdf_4.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(rad_cal_fp_mdf_0_grid,rad_cal_fp_mdf_0_values,'b+');
plot(rad_cal_fp_mc_seq_kdeEvalPositions(1,:),rad_cal_fp_mc_seq_gaussianKdeDensities(1,:),'b-','linewidth',3);
hold
%%plot(rad_val_fp_mdf_0_grid,rad_val_fp_mdf_0_values,'r+');
plot(rad_val_fp_mc_seq_kdeEvalPositions(1,:),rad_val_fp_mc_seq_gaussianKdeDensities(1,:),'r-','linewidth',3);
ylabel('Pdf','fontsize',20);
xlabel('Radiative heat flux [W/m^{2}]','fontsize',20);
title('QoI Pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng talk_rad_cal_val_qoi_mpdf_5.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(rad_cal_fp_cdf_0_grid,100*rad_cal_fp_cdf_0_values,'b-','linewidth',3);
hold
plot(rad_val_fp_cdf_0_grid,100*rad_val_fp_cdf_0_values,'r-','linewidth',3);

a=axis;
%axis([a(1) 1.0 a(3) a(4)]);
ylabel('Cdf (%)','fontsize',20);
xlabel('Radiative heat flux [W/m^{2}]','fontsize',20);
title('QoI Cdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng talk_rad_cal_val_qoi_cdf_6.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot(rad_cal_fp_cdf_0_grid,100*rad_cal_fp_cdf_0_values,'b-','linewidth',3);
%hold
%plot(rad_val_fp_cdf_0_grid,100*rad_val_fp_cdf_0_values,'r-','linewidth',3);

%[m,n] = size(rad_cal_fp_cdf_0_grid);
%epsilon = 0.06;
%plot([0:0.01:1],ones(101,1)*100*(0.+epsilon/2),'k--','linewidth',1);
%plot([0:0.01:1],ones(101,1)*100*(1.-epsilon/2),'k--','linewidth',1);

%a=axis;
%axis([a(1) 1.0 a(3) a(4)]);
%ylabel('Cdf (%)','fontsize',20);
%xlabel('Radiative heat flux [W/m^{2}]','fontsize',20);
%title('Full Picture','fontsize',20);
%grid on;
%set(gca,'fontsize',20);
%legend('Calibration',...
%       'Validation',...
%       'location','northwest');
%print -dpng talk_rad_cal_val_qoi_cdf_0_model_confidence_a.png
%waitforbuttonpress;
%clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot(rad_cal_fp_cdf_0_grid,100*rad_cal_fp_cdf_0_values,'b-','linewidth',3);
%hold
%plot(rad_val_fp_cdf_0_grid,100*rad_val_fp_cdf_0_values,'r-','linewidth',3);

%[m,n] = size(rad_cal_fp_cdf_0_grid);
%epsilon = 0.06;
%plot([0:0.01:1],ones(101,1)*100*(0.+epsilon/2),'k--','linewidth',3);
%plot([0:0.01:1],ones(101,1)*100*(1.-epsilon/2),'k--','linewidth',3);

%a=axis;
%axis([0.15 0.27 2 15]);
%ylabel('Cdf (%)','fontsize',20);
%xlabel('Radiative heat flux [W/m^{2}]','fontsize',20);
%title('Zoom','fontsize',20);
%grid minor;
%set(gca,'fontsize',20);
%legend('Calibration',...
%       'Validation',...
%       'location','northwest');
%print -dpng talk_rad_cal_val_qoi_cdf_1_model_confidence_b.png
%waitforbuttonpress;
%clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot(rad_cal_fp_cdf_0_grid,100*rad_cal_fp_cdf_0_values,'b-','linewidth',3);
%hold
%plot(rad_val_fp_cdf_0_grid,100*rad_val_fp_cdf_0_values,'r-','linewidth',3);

%[m,n] = size(rad_cal_fp_cdf_0_grid);
%probabilityTresholdForFailure = 0.05;
%plot([0:0.01:1],ones(101,1)*100*probabilityTresholdForFailure,'k--','linewidth',1);
%massTresholdForFailure = 0.2;
%plot(ones(101,1)*massTresholdForFailure,[0:1:100],'k--','linewidth',1);

%a=axis;
%axis([a(1) 1.0 a(3) a(4)]);
%ylabel('Cdf (%)','fontsize',20);
%xlabel('Radiative heat flux [W/m^{2}]','fontsize',20);
%title('Full Picture','fontsize',20);
%grid on;
%set(gca,'fontsize',20);
%legend('Calibration',...
%       'Validation',...
%       'location','northwest');
%print -dpng talk_rad_cal_val_qoi_cdf_2_system_confidence_a.png
%waitforbuttonpress;
%clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot(rad_cal_fp_cdf_0_grid,100*rad_cal_fp_cdf_0_values,'b-','linewidth',3);
%hold
%plot(rad_val_fp_cdf_0_grid,100*rad_val_fp_cdf_0_values,'r-','linewidth',3);

%[m,n] = size(rad_cal_fp_cdf_0_grid);
%probabilityTresholdForFailure = 0.05;
%plot([0:0.01:1],ones(101,1)*100*probabilityTresholdForFailure,'k--','linewidth',3);
%massTresholdForFailure = 0.2;
%plot(ones(101,1)*massTresholdForFailure,[0:1:100],'k--','linewidth',3);

%a=axis;
%axis([0.15 0.25 2 8]);
%ylabel('Cdf (%)','fontsize',20);
%xlabel('Radiative heat flux [W/m^{2}]','fontsize',20);
%title('Zoom','fontsize',20);
%grid minor;
%set(gca,'fontsize',20);
%legend('Calibration',...
%       'Validation',...
%       'location','northwest');
%print -dpng talk_rad_cal_val_qoi_cdf_3_system_confidence_b.png
%waitforbuttonpress;
%clf
