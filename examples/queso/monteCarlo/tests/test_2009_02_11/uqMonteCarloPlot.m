cd outputData
mcOutput
mcExtraOutput
cd ..

%fp_mc_seq_corrViaFft_lags = zeros(1,15);
%fp_mc_seq_corrViaFft_initPos0 = zeros(1,15);
%fp_mc_seq_kdeEvalPositions = zeros(1,250);
%fp_mc_seq_gaussianKdeScaleVec = zeros(1);
%fp_mc_seq_gaussianKdeDensities = zeros(1,250);
%fp_mdf_0_grid = zeros(250,1);
%fp_mdf_0_values = zeros(250,1);
%fp_cdf_0_grid = zeros(250,1);
%fp_cdf_0_values = zeros(250,1);

plot(fp_mc_seq_corrViaFft_lags,fp_mc_seq_corrViaFft_initPos0,'-o');
title('QoI autocorrelation');
waitforbuttonpress;
clf;

plot(fp_mc_seq_kdeEvalPositions,fp_mc_seq_gaussianKdeDensities,'-');
title('QoI KDE');
waitforbuttonpress;
clf;

plot(fp_mdf_0_grid,fp_mdf_0_values,'-');
title('QoI mdf');
waitforbuttonpress;
clf;

plot(fp_cdf_0_grid,fp_cdf_0_values,'-');
title('QoI cdf');
waitforbuttonpress;
clf;

