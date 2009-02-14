s1Output
s1ExtraOutput

plot(s1_cal_mdf_0_grid,s1_cal_mdf_0_values,'b+');
hold
plot(s1_cal_mc_filteredChain_kdeEvalPositions(1,:),s1_cal_mc_filteredChain_gaussianKdeDensities(1,:),'bo');
title('Stage 1, pdf of calibrated A: scaled histogram (+) and gaussian kde(o)');
print -dpng s1_verify4_fig1.png
waitforbuttonpress;
clf

plot(s1_cal_mdf_1_grid,s1_cal_mdf_1_values,'b+');
hold
plot(s1_cal_mc_filteredChain_kdeEvalPositions(2,:),s1_cal_mc_filteredChain_gaussianKdeDensities(2,:),'bo');
title('Stage 1, pdf of calibrated E: scaled histogram (+) and gaussian kde(o)');
print -dpng s1_verify4_fig2.png
waitforbuttonpress;
clf

plot(s1_pro_mdf_0_grid,s1_pro_mdf_0_values,'b+');
hold
plot(s1_pro_mc_seq_kdeEvalPositions(1,:),s1_pro_mc_seq_gaussianKdeDensities(1,:),'bo');
title('Stage 1, pdf of QoI t.25: scaled histogram (+) and gaussian kde(o)');
print -dpng s1_verify4_fig3.png
waitforbuttonpress;
clf

plot(s1_pro_cdf_0_grid,s1_pro_cdf_0_values,'bo');
title('Stage 1, cdf of QoI t.25 (via histogram)');
print -dpng s1_verify4_fig4.png
waitforbuttonpress;
clf
