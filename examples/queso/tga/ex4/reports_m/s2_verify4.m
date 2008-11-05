s2Output
s2ExtraOutput

plot(s2_cal_mdf_0_grid,s2_cal_mdf_0_values,'r+');
hold
plot(s2_cal_mc_filteredChain_kdeEvalPositions(1,:),s2_cal_mc_filteredChain_gaussianKdeDensities(1,:),'ro');
title('Stage 2, pdf of calibrated A: scaled histogram (+) and gaussian kde(o)');
print -dpng s2_verify4_fig1.png
waitforbuttonpress;
clf

plot(s2_cal_mdf_1_grid,s2_cal_mdf_1_values,'r+');
hold
plot(s2_cal_mc_filteredChain_kdeEvalPositions(2,:),s2_cal_mc_filteredChain_gaussianKdeDensities(2,:),'ro');
title('Stage 2, pdf of calibrated E: scaled histogram (+) and gaussian kde(o)');
print -dpng s2_verify4_fig2.png
waitforbuttonpress;
clf

plot(s2_pro_mdf_0_grid,s2_pro_mdf_0_values,'r+');
hold
plot(s2_pro_mc_seq_kdeEvalPositions(1,:),s2_pro_mc_seq_gaussianKdeDensities(1,:),'ro');
title('Stage 2, pdf of QoI t.25: scaled histogram (+) and gaussian kde(o)');
print -dpng s2_verify4_fig3.png
waitforbuttonpress;
clf

plot(s2_pro_cdf_0_grid,s2_pro_cdf_0_values,'ro');
title('Stage 2, cdf of QoI t.25 (via histogram)');
print -dpng s2_verify4_fig4.png
waitforbuttonpress;
clf
