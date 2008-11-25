s1Output
s1ExtraOutput
s2Output
s2ExtraOutput

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minPrior0 = 2.40e+11;
maxPrior0 = 2.80e+11;
numPrior0 = 200;
numHorizPts0 = 20;
deltaPrior0 = (maxPrior0 - minPrior0)/numPrior0;
atc_cal_ip_prior_0_grid   = [minPrior0-numHorizPts0*deltaPrior0 : deltaPrior0 : maxPrior0+numHorizPts0*deltaPrior0];
atc_cal_ip_prior_0_values = ones(1,numPrior0+2*numHorizPts0+1)./(maxPrior0-minPrior0);
atc_cal_ip_prior_0_values(1,1:numHorizPts0)     = 0.;
atc_cal_ip_prior_0_values(1,numPrior0+numHorizPts0+2:numPrior0+2*numHorizPts0+1) = 0.;
plot(atc_cal_ip_prior_0_grid,atc_cal_ip_prior_0_values,'b*');
hold

atc_cal_left_vertical_line_x0 = ones(1,11)*minPrior0;
atc_cal_left_vertical_line_y0 = [0 : atc_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : atc_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(atc_cal_left_vertical_line_x0,atc_cal_left_vertical_line_y0,'b--','linewidth',1);

atc_cal_right_vertical_line_x0 = ones(1,11)*maxPrior0;
atc_cal_right_vertical_line_y0 = [0 : atc_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : atc_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(atc_cal_right_vertical_line_x0,atc_cal_right_vertical_line_y0,'b--','linewidth',1);

ylabel('Prior marginal pdf','fontsize',20);
xlabel('A (min^{-1})','fontsize',20);
title('Fig 1: Parameter A: prior marginal pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng talk_atc_cal_prior_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minPrior1 = 1.80e+5;
maxPrior1 = 2.20e+5;
numPrior1 = 200;
numHorizPts1 = 20;
deltaPrior1 = (maxPrior1 - minPrior1)/numPrior1;
atc_cal_ip_prior_1_grid   = [minPrior1-numHorizPts1*deltaPrior1 : deltaPrior1 : maxPrior1+numHorizPts1*deltaPrior1];
atc_cal_ip_prior_1_values = ones(1,numPrior1+2*numHorizPts1+1)./(maxPrior1-minPrior1);
atc_cal_ip_prior_1_values(1,1:numHorizPts1)     = 0.;
atc_cal_ip_prior_1_values(1,numPrior1+numHorizPts1+2:numPrior1+2*numHorizPts1+1) = 0.;
plot(atc_cal_ip_prior_1_grid,atc_cal_ip_prior_1_values,'b*');
hold

atc_cal_left_vertical_line_x1 = ones(1,11)*minPrior1;
atc_cal_left_vertical_line_y1 = [0 : atc_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : atc_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(atc_cal_left_vertical_line_x1,atc_cal_left_vertical_line_y1,'b--','linewidth',1);

atc_cal_right_vertical_line_x1 = ones(1,11)*maxPrior1;
atc_cal_right_vertical_line_y1 = [0 : atc_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : atc_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(atc_cal_right_vertical_line_x1,atc_cal_right_vertical_line_y1,'b--','linewidth',1);

ylabel('Prior marginal pdf','fontsize',20);
xlabel('E (J/mol)','fontsize',20);
title('Fig 2: Parameter E: prior marginal pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng talk_atc_cal_prior_1.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(atc_cal_ip_mdf_0_grid,atc_cal_ip_mdf_0_values,'b+');
plot(atc_cal_ip_mc_filteredChain_kdeEvalPositions(1,:),atc_cal_ip_mc_filteredChain_gaussianKdeDensities(1,:),'b-','linewidth',3);
hold
%%plot(atc_val_ip_mdf_0_grid,atc_val_ip_mdf_0_values,'r+');
plot(atc_val_ip_mc_filteredChain_kdeEvalPositions(1,:),atc_val_ip_mc_filteredChain_gaussianKdeDensities(1,:),'r-','linewidth',3);
ylabel('Posterior marginal pdf','fontsize',20);
xlabel('A (min^{-1})','fontsize',20);
title('Fig 3: A: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','southwest');

plot(atc_cal_ip_prior_0_grid,atc_cal_ip_prior_0_values,'b*');

atc_cal_left_vertical_line_x0 = ones(1,11)*minPrior0;
atc_cal_left_vertical_line_y0 = [0 : atc_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : atc_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(atc_cal_left_vertical_line_x0,atc_cal_left_vertical_line_y0,'b--','linewidth',1);

atc_cal_right_vertical_line_x0 = ones(1,11)*maxPrior0;
atc_cal_right_vertical_line_y0 = [0 : atc_cal_ip_prior_0_values(1,numHorizPts0+1)/10 : atc_cal_ip_prior_0_values(1,numHorizPts0+1)];
plot(atc_cal_right_vertical_line_x0,atc_cal_right_vertical_line_y0,'b--','linewidth',1);

print -dpng talk_atc_cal_val_post_mpdf_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(atc_cal_ip_mdf_1_grid,atc_cal_ip_mdf_1_values,'b+');
plot(atc_cal_ip_mc_filteredChain_kdeEvalPositions(2,:),atc_cal_ip_mc_filteredChain_gaussianKdeDensities(2,:),'b-','linewidth',3);
hold
%%plot(atc_val_ip_mdf_1_grid,atc_val_ip_mdf_1_values,'r+');
plot(atc_val_ip_mc_filteredChain_kdeEvalPositions(2,:),atc_val_ip_mc_filteredChain_gaussianKdeDensities(2,:),'r-','linewidth',3);
ylabel('Pdf','fontsize',20);
xlabel('E (J/mol)','fontsize',20);
title('Fig 4: E: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');

plot(atc_cal_ip_prior_1_grid,atc_cal_ip_prior_1_values,'b*');

atc_cal_left_vertical_line_x1 = ones(1,11)*minPrior1;
atc_cal_left_vertical_line_y1 = [0 : atc_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : atc_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(atc_cal_left_vertical_line_x1,atc_cal_left_vertical_line_y1,'b--','linewidth',1);

atc_cal_right_vertical_line_x1 = ones(1,11)*maxPrior1;
atc_cal_right_vertical_line_y1 = [0 : atc_cal_ip_prior_1_values(1,numHorizPts1+1)/10 : atc_cal_ip_prior_1_values(1,numHorizPts1+1)];
plot(atc_cal_right_vertical_line_x1,atc_cal_right_vertical_line_y1,'b--','linewidth',1);

print -dpng talk_atc_cal_val_post_mpdf_1.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(atc_cal_fp_mdf_0_grid,atc_cal_fp_mdf_0_values,'b+');
plot(atc_cal_fp_mc_seq_kdeEvalPositions(1,:),atc_cal_fp_mc_seq_gaussianKdeDensities(1,:),'b-','linewidth',3);
hold
%%plot(atc_val_fp_mdf_0_grid,atc_val_fp_mdf_0_values,'r+');
plot(atc_val_fp_mc_seq_kdeEvalPositions(1,:),atc_val_fp_mc_seq_gaussianKdeDensities(1,:),'r-','linewidth',3);
ylabel('Pdf','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('Fig 5: QoI Pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng talk_atc_cal_val_qoi_mpdf_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(atc_cal_fp_cdf_0_grid,100*atc_cal_fp_cdf_0_values,'b-','linewidth',3);
hold
plot(atc_val_fp_cdf_0_grid,100*atc_val_fp_cdf_0_values,'r-','linewidth',3);

a=axis;
%axis([a(1) 1.0 a(3) a(4)]);
ylabel('Cdf (%)','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('Fig 6: QoI Cdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng talk_atc_cal_val_qoi_cdf_0.png
waitforbuttonpress;
clf
