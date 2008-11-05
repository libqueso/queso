stage0Output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots on parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 1;
if a > 0

plot(val_stage_0_cal_bmcdc_0_chain_centersOfHistBins(1,:),val_stage_0_cal_bmcdc_0_chain_histBins(1,:),'o');
waitforbuttonpress; %sleep(3);
clf;

plot(val_stage_0_cal_bmcdc_0_chain_centersOfHistBins(2,:),val_stage_0_cal_bmcdc_0_chain_histBins(2,:),'o');
waitforbuttonpress;
clf;

plot(val_stage_0_cal_bmcdc_0_chain_kdeEvalPositions(1,:),val_stage_0_cal_bmcdc_0_chain_gaussianKdeDensities(1,:),'o'); 
waitforbuttonpress;
clf;

plot(val_stage_0_cal_bmcdc_0_chain_kdeEvalPositions(2,:),val_stage_0_cal_bmcdc_0_chain_gaussianKdeDensities(2,:),'o');
waitforbuttonpress;
clf;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots on qois
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(val_stage_0_pro_mcdc_seq_centersOfHistBins(1,:),val_stage_0_pro_mcdc_seq_histBins(1,:),'o');
waitforbuttonpress;
clf;

plot(val_stage_0_pro_mcdc_seq_kdeEvalPositions(1,:),val_stage_0_pro_mcdc_seq_gaussianKdeDensities(1,:),'o'); 


