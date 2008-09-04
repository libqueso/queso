calib1
[Pxx,W]=pwelch(queso_calib_0_chain);

plot(1:2049,queso_calib_0_chain_psd_initPos0,'bo');
hold
plot(1:2049,Pxx,'r*');
title('Difference between components of uq_psd and matlab_pwelch');

