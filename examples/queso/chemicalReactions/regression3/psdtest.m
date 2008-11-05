calib1

wind=hamming(3640,'symmetric');
[Pxx3,W]=pwelch(queso_calib_0_chain,wind,1820,4096,'twosided');

plot(1,queso_calib_0_chain_psd_initPos0(1),'bo');
hold
plot(1,Pxx3(1),'r*');
legend('QUESO-MCMC Tool',...
       'Matlab pwelch',...
       'Location','NorthEast');
plot(2:2049,queso_calib_0_chain_psd_initPos0(2:2049),'bo');
plot(2:2049,Pxx3(2:2049),'r*');
title('Difference between components of uq_psd and matlab_pwelch');

%waitforbuttonpress;
%clf
%plot(2:2048,queso_calib_0_chain_psd_initPos0(1,2:2048)./Pxx3(2:2048,1)','bo'); %'

waitforbuttonpress;
clf
plot(1:2049,queso_calib_0_chain_psd_initPos0(1,1:2049)./Pxx3(1:2049,1)','bo'); %'
title('Ratio between components of uq_psd and matlab_pwelch');

%[Pxx,W]=pwelch(queso_calib_0_chain);

%[Pxx2,W]=pwelch(queso_calib_0_chain,3640,1820,4096,'onesided');
%plot(1:2049,Pxx-Pxx2,'r*');

%wind=hamming(3640,'symmetric');
%[Pxx2,W]=pwelch(queso_calib_0_chain,wind,1820,4096,'onesided');
%plot(1:2049,Pxx-Pxx2,'r*');

