s1Output
s1ExtraOutput
s2Output
s2ExtraOutput

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minPrior0 = 2.40e+11;
maxPrior0 = 2.80e+11;
numPrior0 = 200;
deltaPrior0 = (maxPrior0 - minPrior0)/numPrior0;
s1_cal_prior_0_grid   = [minPrior0-5*deltaPrior0 : deltaPrior0 : maxPrior0+5*deltaPrior0];
s1_cal_prior_0_values = ones(1,numPrior0+11)./(maxPrior0-minPrior0);
s1_cal_prior_0_values(1,1:5)     = 0.;
s1_cal_prior_0_values(1,numPrior0+7:numPrior0+11) = 0.;
plot(s1_cal_prior_0_grid,s1_cal_prior_0_values,'b*');
ylabel('Prior pdf');
xlabel('A');
title('Stage 1, parameter A: prior pdf');
print -dpng s1_prior_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minPrior1 = 1.80e+5;
maxPrior1 = 2.20e+5;
numPrior1 = 200;
deltaPrior1 = (maxPrior1 - minPrior1)/numPrior1;
s1_cal_prior_1_grid   = [minPrior1-5*deltaPrior1 : deltaPrior1 : maxPrior1+5*deltaPrior1];
s1_cal_prior_1_values = ones(1,numPrior1+11)./(maxPrior1-minPrior1);
s1_cal_prior_1_values(1,1:5)     = 0.;
s1_cal_prior_1_values(1,numPrior1+7:numPrior1+11) = 0.;
plot(s1_cal_prior_1_grid,s1_cal_prior_1_values,'b*');
ylabel('Prior pdf');
xlabel('E');
title('Stage 1, parameter E: prior pdf');
print -dpng s1_prior_1.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(s1_cal_mc_chain_corrViaFft_lags,s1_cal_mc_chain_corrViaFft_initPos0(1,:),'b-o');
hold
plot(s2_cal_mc_chain_corrViaFft_lags,s2_cal_mc_chain_corrViaFft_initPos0(1,:),'r-o');
a=axis;
axis([a(1) a(2) 0. a(4)]);
ylabel('Autocorrelation');
xlabel('Lag');
title('Stages 1 (blue) and 2 (red), parameter A: autocorrelation of "raw" Markov Chain');
print -dpng s1_s2_post_corr_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(s1_cal_mc_chain_corrViaFft_lags,s1_cal_mc_chain_corrViaFft_initPos0(2,:),'b-o');
hold
plot(s2_cal_mc_chain_corrViaFft_lags,s2_cal_mc_chain_corrViaFft_initPos0(2,:),'r-o');
a=axis;
axis([a(1) a(2) 0. a(4)]);
ylabel('Autocorrelation');
xlabel('Lag');
title('Stages 1 (blue) and 2 (red), parameter E: autocorrelation of "raw" Markov Chain');
print -dpng s1_s2_post_corr_1.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(s1_cal_prior_0_grid,s1_cal_prior_0_values,'b*');
hold
plot(s1_cal_mdf_0_grid,s1_cal_mdf_0_values,'b+');
plot(s1_cal_mc_filteredChain_kdeEvalPositions(1,:),s1_cal_mc_filteredChain_gaussianKdeDensities(1,:),'bo');
plot(s2_cal_mdf_0_grid,s2_cal_mdf_0_values,'r+');
plot(s2_cal_mc_filteredChain_kdeEvalPositions(1,:),s2_cal_mc_filteredChain_gaussianKdeDensities(1,:),'ro');
ylabel('Pdf');
xlabel('A');
title('Stages 1 (blue) and 2 (red), param A: prior (*), post histogram (+) and post gaussian kde(o)');
print -dpng s1_s2_post_mdf_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(s1_cal_prior_1_grid,s1_cal_prior_1_values,'b*');
hold
plot(s1_cal_mdf_1_grid,s1_cal_mdf_1_values,'b+');
plot(s1_cal_mc_filteredChain_kdeEvalPositions(2,:),s1_cal_mc_filteredChain_gaussianKdeDensities(2,:),'bo');
plot(s2_cal_mdf_1_grid,s2_cal_mdf_1_values,'r+');
plot(s2_cal_mc_filteredChain_kdeEvalPositions(2,:),s2_cal_mc_filteredChain_gaussianKdeDensities(2,:),'ro');
ylabel('Pdf');
xlabel('E');
title('Stages 1 (blue) and 2 (red), param E: prior (*), post histogram (+) and post gaussian kde(o)');
print -dpng s1_s2_post_mdf_1.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(s1_pro_mc_seq_corrViaFft_lags,s1_pro_mc_seq_corrViaFft_initPos0(1,:),'b-o');
hold
plot(s2_pro_mc_seq_corrViaFft_lags,s2_pro_mc_seq_corrViaFft_initPos0(1,:),'r-o');
%a=axis;
%axis([a(1) a(2) 0. a(4)]);
ylabel('Autocorrelation');
xlabel('Lag');
title('Stages 1 (blue) and 2 (red), qoi: autocorrelation of Monte Carlo output "chain"');
print -dpng s1_s2_qoi_corr_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(s1_pro_mdf_0_grid,s1_pro_mdf_0_values,'b+');
hold
plot(s1_pro_mc_seq_kdeEvalPositions(1,:),s1_pro_mc_seq_gaussianKdeDensities(1,:),'bo');
plot(s2_pro_mdf_0_grid,s2_pro_mdf_0_values,'r+');
plot(s2_pro_mc_seq_kdeEvalPositions(1,:),s2_pro_mc_seq_gaussianKdeDensities(1,:),'ro');
ylabel('Pdf');
xlabel('QoI (seconds)');
title('Stages 1 (blue) and 2 (red), qoi: scaled histogram (+) and gaussian kde(o)');
print -dpng s1_s2_qoi_mdf_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(s1_pro_cdf_0_grid,s1_pro_cdf_0_values,'bo');
hold
plot(s2_pro_cdf_0_grid,s2_pro_cdf_0_values,'ro');
ylabel('Cdf');
xlabel('QoI (seconds)');
title('Stages 1 (blue) and 2 (red), qoi: cdf via histogram');
print -dpng s1_s2_qoi_cdf_0.png
waitforbuttonpress;
clf

