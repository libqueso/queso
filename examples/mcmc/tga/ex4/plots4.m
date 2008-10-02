s1Output
s2Output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stage I: Plots on correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(s1_cal_mc_chain_corrViaFft_lags(1,:),s1_cal_mc_chain_corrViaFft_initPos0(1,:),'bo-');
title('Stage I: Autocorrelation plot for parameter A, on original chain of 1,048,576 positions');
print -dpng s1_corrA.png
waitforbuttonpress; %sleep(3);
clf;

%plot(s1_cal_mc_chain_corrViaFft_lags(2,:),s1_cal_mc_chain_corrViaFft_initPos0(2,:),'bo-');
%title('Stage I: Autocorrelation plot for parameter E, on original chain of 1,048,576 positions');
%print -dpng s1_corrE.png
%waitforbuttonpress; %sleep(3);
%clf;

plot(s1_pro_mc_seq_corrViaFft_lags(1,:),s1_pro_mc_seq_corrViaFft_initPos0(1,:),'bo-');
title('Stage I: Autocorrelation plot for qoi, on sequence of 52,429 samples');
print -dpng s1_corrQ.png
waitforbuttonpress; %sleep(3);
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stage I: Plots on parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 1;
if a > 0

plot(s1_cal_mc_filteredChain_centersOfHistBins(1,:),s1_cal_mc_filteredChain_histBins(1,:),'bo');
title('Stage I: Posterior histogram for parameter A');
print -dpng s1_histA.png
waitforbuttonpress; %sleep(3);
clf;

plot(s1_cal_mc_filteredChain_centersOfHistBins(2,:),s1_cal_mc_filteredChain_histBins(2,:),'bo');
title('Stage I: Posterior histogram for parameter E');
print -dpng s1_histE.png
waitforbuttonpress;
clf;

plot(s1_cal_mc_filteredChain_kdeEvalPositions(1,:),s1_cal_mc_filteredChain_gaussianKdeDensities(1,:),'bo'); 
title('Stage I: Posterior density for parameter A');
print -dpng s1_densA.png
waitforbuttonpress;
clf;

plot(s1_cal_mc_filteredChain_kdeEvalPositions(2,:),s1_cal_mc_filteredChain_gaussianKdeDensities(2,:),'bo');
title('Stage I: Posterior density for parameter E');
print -dpng s1_densE.png
waitforbuttonpress;
clf;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stage I: Plots on qois
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(s1_pro_mc_seq_centersOfHistBins(1,:),s1_pro_mc_seq_histBins(1,:),'bo');
title('Stage I: Histogram for qoi = time for .25 residual mass');
print -dpng s1_histQ.png
waitforbuttonpress;
clf;

plot(s1_pro_mc_seq_kdeEvalPositions(1,:),s1_pro_mc_seq_gaussianKdeDensities(1,:),'bo'); 
hold
plot(s2_pro_mc_seq_kdeEvalPositions(1,:),s2_pro_mc_seq_gaussianKdeDensities(1,:),'ro'); 
title('Stages I (blues) and II (red): Density for qoi = time for .25 residual mass');
print -dpng s1_densQ.png

