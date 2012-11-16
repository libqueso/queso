cd outputData
sfpOutput_sub0
%sfpExtraOutput_sub0
%cd ..

plot(fp_mc_QoiSeq_corrViaFftLags_sub0,fp_mc_QoiSeq_corrViaFftInitPos0_sub0,'-o');
title('QoI autocorrelation');
print -dpng fig1.png
waitforbuttonpress;
clf;

plot(fp_mc_QoiSeq_unifGkdePosits_sub0,fp_mc_QoiSeq_unifGkdeValues_sub0,'-');
title('QoI KDE');
print -dpng fig2.png
waitforbuttonpress;
clf;

plot(fp_unifQoiCdf_0_grid_sub0,fp_unifQoiCdf_0_values_sub0,'-');
title('QoI cdf');
print -dpng fig3.png
waitforbuttonpress;
clf;

