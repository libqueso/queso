s1Output
s1ExtraOutput
s2Output
s2ExtraOutput

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%plot(s1_pro_mdf_0_grid,s1_pro_mdf_0_values,'b+');
plot(s1_pro_mc_seq_kdeEvalPositions(1,:),s1_pro_mc_seq_gaussianKdeDensities(1,:),'b-');
hold
%%plot(s2_pro_mdf_0_grid,s2_pro_mdf_0_values,'r+');
plot(s2_pro_mc_seq_kdeEvalPositions(1,:),s2_pro_mc_seq_gaussianKdeDensities(1,:),'r-');
ylabel('Pdf');
xlabel('QoI (seconds)');
title('Stages 1 (blue) and 2 (red), qoi: gaussian kdes(o)');
print -dpng s1_s2_qoi_mdf_0_250Kmin.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(s1_pro_cdf_0_grid,s1_pro_cdf_0_values,'b-');
hold
plot(s2_pro_cdf_0_grid,s2_pro_cdf_0_values,'r-');

[m,n] = size(s1_pro_cdf_0_grid);
epsilon = 0.08;
plot(s1_pro_cdf_0_grid,ones(m,n)*(0.+epsilon/2),'k--');
plot(s1_pro_cdf_0_grid,ones(m,n)*(1.-epsilon/2),'k--');

ylabel('Cdf');
xlabel('QoI (seconds)');
title('Stages 1 (blue) and 2 (red), qoi: cdf (via histogram) and horizontal lines at \epsilon/2 and 1-\epsilon/2 for \epsilon=0.08');
print -dpng s1_s2_qoi_cdf_0_250Kmin.png
waitforbuttonpress;
clf

