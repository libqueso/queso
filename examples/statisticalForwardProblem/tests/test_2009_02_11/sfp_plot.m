cd outputData
sfpOutput_sub0
%sfpExtraOutput_sub0
%cd ..

plot(fp_mc_QoiSeq_corrViaFftLags_sub0,fp_mc_QoiSeq_corrViaFftInitPos0_sub0,'-o');
title('QoI autocorrelation');
waitforbuttonpress;
clf;

plot(fp_mc_QoiSeq_GkdePosits_sub0,fp_mc_QoiSeq_GkdeValues_sub0,'-');
title('QoI KDE');
waitforbuttonpress;
clf;

plot(fp_QoiMdf_0_grid_sub0,fp_QoiMdf_0_values_sub0,'-');
title('QoI mdf');
waitforbuttonpress;
clf;

plot(fp_unifQoiCdf_0_grid_sub0,fp_unifQoiCdf_0_values_sub0,'-');
title('QoI cdf');
waitforbuttonpress;
clf;

