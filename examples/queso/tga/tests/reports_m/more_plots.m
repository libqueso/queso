cd results_5_25
s1Output
s1ExtraOutput
s2Output
s2ExtraOutput

orig_s1_pro_mc_seq_kdeEvalPositions     = s1_pro_mc_seq_kdeEvalPositions;
orig_s1_pro_mc_seq_gaussianKdeDensities = s1_pro_mc_seq_gaussianKdeDensities;
orig_s2_pro_mc_seq_kdeEvalPositions     = s2_pro_mc_seq_kdeEvalPositions;
orig_s2_pro_mc_seq_gaussianKdeDensities = s2_pro_mc_seq_gaussianKdeDensities;

orig_s1_pro_cdf_0_grid   = s1_pro_cdf_0_grid;
orig_s1_pro_cdf_0_values = s1_pro_cdf_0_values;
orig_s2_pro_cdf_0_grid   = s2_pro_cdf_0_grid;
orig_s2_pro_cdf_0_values = s2_pro_cdf_0_values;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ../results_5_25_50
s1Output
s1ExtraOutput
s2Output
s2ExtraOutput

more_s1_pro_mc_seq_kdeEvalPositions     = s1_pro_mc_seq_kdeEvalPositions;
more_s1_pro_mc_seq_gaussianKdeDensities = s1_pro_mc_seq_gaussianKdeDensities;
more_s2_pro_mc_seq_kdeEvalPositions     = s2_pro_mc_seq_kdeEvalPositions;
more_s2_pro_mc_seq_gaussianKdeDensities = s2_pro_mc_seq_gaussianKdeDensities;

more_s1_pro_cdf_0_grid   = s1_pro_cdf_0_grid;
more_s1_pro_cdf_0_values = s1_pro_cdf_0_values;
more_s2_pro_cdf_0_grid   = s2_pro_cdf_0_grid;
more_s2_pro_cdf_0_values = s2_pro_cdf_0_values;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(orig_s1_pro_mc_seq_kdeEvalPositions(1,:),orig_s1_pro_mc_seq_gaussianKdeDensities(1,:),'b--');
hold
plot(orig_s2_pro_mc_seq_kdeEvalPositions(1,:),orig_s2_pro_mc_seq_gaussianKdeDensities(1,:),'r--');
plot(more_s1_pro_mc_seq_kdeEvalPositions(1,:),more_s1_pro_mc_seq_gaussianKdeDensities(1,:),'b-');
plot(more_s2_pro_mc_seq_kdeEvalPositions(1,:),more_s2_pro_mc_seq_gaussianKdeDensities(1,:),'r-');
ylabel('Pdf','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('QoI: pdf (calibration \equiv blue; validation \equiv red)','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Cal. 5, 25 K/min',...
       'Val. 100 K/min',...
       'Cal. 5, 25, 50 K/min',...
       'Val. 100 K/min',...
       'location','northeast');
print -dpng more_s1_s2_qoi_mpdf_0.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(orig_s1_pro_cdf_0_grid,100*orig_s1_pro_cdf_0_values,'b--');
hold
plot(orig_s2_pro_cdf_0_grid,100*orig_s2_pro_cdf_0_values,'r--');
plot(more_s1_pro_cdf_0_grid,100*more_s1_pro_cdf_0_values,'b-');
plot(more_s2_pro_cdf_0_grid,100*more_s2_pro_cdf_0_values,'r-');

[m,n] = size(orig_s1_pro_cdf_0_grid);
probabilityTresholdForFailure = 0.05;
plot([0:0.01:1],ones(101,1)*100*probabilityTresholdForFailure,'k--','linewidth',2);
massTresholdForFailure = 0.2;
plot(ones(101,1)*massTresholdForFailure,[0:1:100],'k--','linewidth',2);

a=axis;
axis([a(1) 0.6 a(3) a(4)]);
ylabel('Cdf (%)','fontsize',20);
xlabel('Mass fraction remaining at t=3.9s','fontsize',20);
title('QoI: cdf (calibration \equiv blue; validation \equiv red)','fontsize',20);
grid on;
set(gca,'fontsize',20);
legend('Cal. 5, 25 K/min',...
       'Val. 100 K/min',...
       'Cal. 5, 25, 50 K/min',...
       'Val. 100 K/min',...
       'location','northwest');
print -dpng more_s1_s2_qoi_cdf_0.png
waitforbuttonpress;
clf

