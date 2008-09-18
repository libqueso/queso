stage0Output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots on parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 1;
if a > 0

plot(val_stage_0_cal_bmc_0_chain_centersOfHistBins(1,:),val_stage_0_cal_bmc_0_chain_histBins(1,:),'o');
sleep(3); %waitforbuttonpress; %sleep(3);
clf;

plot(val_stage_0_cal_bmc_0_chain_centersOfHistBins(2,:),val_stage_0_cal_bmc_0_chain_histBins(2,:),'o');
sleep(3); %waitforbuttonpress;
clf;

plot(val_stage_0_cal_bmc_0_chain_kdeEvalPositions(1,:),val_stage_0_cal_bmc_0_chain_gaussianKdeDensities(1,:),'o'); 
sleep(3); %waitforbuttonpress;
clf;

plot(val_stage_0_cal_bmc_0_chain_kdeEvalPositions(2,:),val_stage_0_cal_bmc_0_chain_gaussianKdeDensities(2,:),'o');
sleep(3); %waitforbuttonpress;
clf;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots on qois
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(val_stage_0_pro_mc_seq_centersOfHistBins(1,:),val_stage_0_pro_mc_seq_histBins(1,:),'o');
sleep(3); %waitforbuttonpress;
clf;

plot(val_stage_0_pro_mc_seq_kdeEvalPositions(1,:),val_stage_0_pro_mc_seq_gaussianKdeDensities(1,:),'o'); 


