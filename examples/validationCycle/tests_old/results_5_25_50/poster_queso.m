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
s1_cal_prior_0_grid   = [minPrior0-numHorizPts0*deltaPrior0 : deltaPrior0 : maxPrior0+numHorizPts0*deltaPrior0];
s1_cal_prior_0_values = ones(1,numPrior0+2*numHorizPts0+1)./(maxPrior0-minPrior0);
s1_cal_prior_0_values(1,1:numHorizPts0)     = 0.;
s1_cal_prior_0_values(1,numPrior0+numHorizPts0+2:numPrior0+2*numHorizPts0+1) = 0.;
plot(s1_cal_prior_0_grid,s1_cal_prior_0_values,'b*');
hold

s1_left_vertical_line_x0 = ones(1,11)*minPrior0;
s1_left_vertical_line_y0 = [0 : s1_cal_prior_0_values(1,numHorizPts0+1)/10 : s1_cal_prior_0_values(1,numHorizPts0+1)];
plot(s1_left_vertical_line_x0,s1_left_vertical_line_y0,'b--','linewidth',1);

s1_right_vertical_line_x0 = ones(1,11)*maxPrior0;
s1_right_vertical_line_y0 = [0 : s1_cal_prior_0_values(1,numHorizPts0+1)/10 : s1_cal_prior_0_values(1,numHorizPts0+1)];
plot(s1_right_vertical_line_x0,s1_right_vertical_line_y0,'b--','linewidth',1);

ylabel('Prior marginal pdf','fontsize',20);
xlabel('A (min^{-1})','fontsize',20);
title('Parameter A: prior marginal pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng poster_queso_s1_prior_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minPrior1 = 1.80e+5;
maxPrior1 = 2.20e+5;
numPrior1 = 200;
numHorizPts1 = 20;
deltaPrior1 = (maxPrior1 - minPrior1)/numPrior1;
s1_cal_prior_1_grid   = [minPrior1-numHorizPts1*deltaPrior1 : deltaPrior1 : maxPrior1+numHorizPts1*deltaPrior1];
s1_cal_prior_1_values = ones(1,numPrior1+2*numHorizPts1+1)./(maxPrior1-minPrior1);
s1_cal_prior_1_values(1,1:numHorizPts1)     = 0.;
s1_cal_prior_1_values(1,numPrior1+numHorizPts1+2:numPrior1+2*numHorizPts1+1) = 0.;
plot(s1_cal_prior_1_grid,s1_cal_prior_1_values,'b*');
hold

s1_left_vertical_line_x1 = ones(1,11)*minPrior1;
s1_left_vertical_line_y1 = [0 : s1_cal_prior_1_values(1,numHorizPts1+1)/10 : s1_cal_prior_1_values(1,numHorizPts1+1)];
plot(s1_left_vertical_line_x1,s1_left_vertical_line_y1,'b--','linewidth',1);

s1_right_vertical_line_x1 = ones(1,11)*maxPrior1;
s1_right_vertical_line_y1 = [0 : s1_cal_prior_1_values(1,numHorizPts1+1)/10 : s1_cal_prior_1_values(1,numHorizPts1+1)];
plot(s1_right_vertical_line_x1,s1_right_vertical_line_y1,'b--','linewidth',1);

ylabel('Prior marginal pdf','fontsize',20);
xlabel('E (J/mol)','fontsize',20);
title('Parameter E: prior marginal pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
print -dpng poster_queso_s1_prior_1.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(s1_cal_mdf_0_grid,s1_cal_mdf_0_values,'b+');
plot(s1_cal_mc_filteredChain_kdeEvalPositions(1,:),s1_cal_mc_filteredChain_gaussianKdeDensities(1,:),'b-','linewidth',3);
hold
%%plot(s2_cal_mdf_0_grid,s2_cal_mdf_0_values,'r+');
plot(s2_cal_mc_filteredChain_kdeEvalPositions(1,:),s2_cal_mc_filteredChain_gaussianKdeDensities(1,:),'r-','linewidth',3);
ylabel('Posterior marginal pdf','fontsize',20);
xlabel('A (min^{-1})','fontsize',20);
title('A: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','southwest');

plot(s1_cal_prior_0_grid,s1_cal_prior_0_values,'b*');

s1_left_vertical_line_x0 = ones(1,11)*minPrior0;
s1_left_vertical_line_y0 = [0 : s1_cal_prior_0_values(1,numHorizPts0+1)/10 : s1_cal_prior_0_values(1,numHorizPts0+1)];
plot(s1_left_vertical_line_x0,s1_left_vertical_line_y0,'b--','linewidth',1);

s1_right_vertical_line_x0 = ones(1,11)*maxPrior0;
s1_right_vertical_line_y0 = [0 : s1_cal_prior_0_values(1,numHorizPts0+1)/10 : s1_cal_prior_0_values(1,numHorizPts0+1)];
plot(s1_right_vertical_line_x0,s1_right_vertical_line_y0,'b--','linewidth',1);

print -dpng poster_queso_s1_s2_post_mpdf_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(s1_cal_mdf_1_grid,s1_cal_mdf_1_values,'b+');
plot(s1_cal_mc_filteredChain_kdeEvalPositions(2,:),s1_cal_mc_filteredChain_gaussianKdeDensities(2,:),'b-','linewidth',3);
hold
%%plot(s2_cal_mdf_1_grid,s2_cal_mdf_1_values,'r+');
plot(s2_cal_mc_filteredChain_kdeEvalPositions(2,:),s2_cal_mc_filteredChain_gaussianKdeDensities(2,:),'r-','linewidth',3);
ylabel('Pdf','fontsize',20);
xlabel('E (J/mol)','fontsize',20);
title('E: prior(*) and posterior (-) marginals','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');

plot(s1_cal_prior_1_grid,s1_cal_prior_1_values,'b*');

s1_left_vertical_line_x1 = ones(1,11)*minPrior1;
s1_left_vertical_line_y1 = [0 : s1_cal_prior_1_values(1,numHorizPts1+1)/10 : s1_cal_prior_1_values(1,numHorizPts1+1)];
plot(s1_left_vertical_line_x1,s1_left_vertical_line_y1,'b--','linewidth',1);

s1_right_vertical_line_x1 = ones(1,11)*maxPrior1;
s1_right_vertical_line_y1 = [0 : s1_cal_prior_1_values(1,numHorizPts1+1)/10 : s1_cal_prior_1_values(1,numHorizPts1+1)];
plot(s1_right_vertical_line_x1,s1_right_vertical_line_y1,'b--','linewidth',1);

print -dpng poster_queso_s1_s2_post_mpdf_1.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(s1_pro_mdf_0_grid,s1_pro_mdf_0_values,'b+');
plot(s1_pro_mc_seq_kdeEvalPositions(1,:),s1_pro_mc_seq_gaussianKdeDensities(1,:),'b-','linewidth',3);
hold
%%plot(s2_pro_mdf_0_grid,s2_pro_mdf_0_values,'r+');
plot(s2_pro_mc_seq_kdeEvalPositions(1,:),s2_pro_mc_seq_gaussianKdeDensities(1,:),'r-','linewidth',3);
ylabel('Pdf','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
%%title('QoI Pdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng poster_queso_s1_s2_qoi_mpdf_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(s1_pro_cdf_0_grid,100*s1_pro_cdf_0_values,'b-','linewidth',3);
hold
plot(s2_pro_cdf_0_grid,100*s2_pro_cdf_0_values,'r-','linewidth',3);

[m,n] = size(s1_pro_cdf_0_grid);
epsilon = 0.06;
plot([0:0.01:1],ones(101,1)*100*(0.+epsilon/2),'k--','linewidth',2);
plot([0:0.01:1],ones(101,1)*100*(1.-epsilon/2),'k--','linewidth',2);

[m,n] = size(s1_pro_cdf_0_grid);
%probabilityTresholdForFailure = 0.05;
%plot([0:0.01:1],ones(101,1)*100*probabilityTresholdForFailure,'k--','linewidth',2);
massTresholdForFailure = 0.2;
plot(ones(101,1)*massTresholdForFailure,[0:1:100],'k--','linewidth',2);

a=axis;
axis([0.1 0.6 a(3) a(4)]);
ylabel('Cdf (%)','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
%%title('QoI Cdf','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Calibration',...
       'Validation',...
       'location','northwest');
print -dpng poster_queso_s1_s2_qoi_cdf_0.png
waitforbuttonpress;
clf


