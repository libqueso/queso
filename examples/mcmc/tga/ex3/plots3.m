cpOutput

plot(cal_mc_chain_corrViaFft_lags(1,:),cal_mc_chain_corrViaFft_initPos0(1,:),'o-');
title('Autocorrelation plot for parameter A, on original chain of 1,048,576 positions');
print -dpng corrA.png
waitforbuttonpress; %sleep(3);
clf;

%plot(cal_mc_chain_corrViaFft_lags(2,:),cal_mc_chain_corrViaFft_initPos0(2,:),'o-');
%title('Autocorrelation plot for parameter E, on original chain of 1,048,576 positions');
%print -dpng corrE.png
%waitforbuttonpress; %sleep(3);
%clf;

plot(pro_mc_seq_corrViaFft_lags(1,:),pro_mc_seq_corrViaFft_initPos0(1,:),'o-');
title('Autocorrelation plot for qoi, on sequence of 52,429 samples');
print -dpng corrQ.png
waitforbuttonpress; %sleep(3);
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots on parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 1;
if a > 0

plot(cal_mc_filteredChain_centersOfHistBins(1,:),cal_mc_filteredChain_histBins(1,:),'o');
title('Posterior histogram for parameter A');
print -dpng histA.png
waitforbuttonpress; %sleep(3);
clf;

plot(cal_mc_filteredChain_centersOfHistBins(2,:),cal_mc_filteredChain_histBins(2,:),'o');
title('Posterior histogram for parameter E');
print -dpng histE.png
waitforbuttonpress;
clf;

plot(cal_mc_filteredChain_kdeEvalPositions(1,:),cal_mc_filteredChain_gaussianKdeDensities(1,:),'o'); 
title('Posterior density for parameter A');
print -dpng densA.png
waitforbuttonpress;
clf;

plot(cal_mc_filteredChain_kdeEvalPositions(2,:),cal_mc_filteredChain_gaussianKdeDensities(2,:),'o');
title('Posterior density for parameter E');
print -dpng densE.png
waitforbuttonpress;
clf;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots on qois
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(pro_mc_seq_centersOfHistBins(1,:),pro_mc_seq_histBins(1,:),'o');
title('Histogram for qoi = time for .25 residual mass');
print -dpng histQ.png
waitforbuttonpress;
clf;

plot(pro_mc_seq_kdeEvalPositions(1,:),pro_mc_seq_gaussianKdeDensities(1,:),'o'); 
title('Density for qoi = time for .25 residual mass');
print -dpng densQ.png


