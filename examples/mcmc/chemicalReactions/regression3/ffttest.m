calib1
abc=fft(queso_calib_0_chain,2048);

plot(1:2048,real(queso_calib_0_chain_fft_initPos0),'bo');
hold
plot(1:2048,real(abc)','r*'); %'
plot(1:2048,real(queso_calib_0_chain_fft_initPos0)-real(abc)','go'); %'
title('Difference between real components of uq_fft and matlab_fft');

waitforbuttonpress;
clf;
plot(1:2048,imag(queso_calib_0_chain_fft_initPos0),'bo');
hold
plot(1:2048,imag(abc)','r*'); %'
plot(1:2048,imag(queso_calib_0_chain_fft_initPos0)-imag(abc)','go'); %'
title('Difference between imag components of uq_fft and matlab_fft');

waitforbuttonpress;
clf;
plot(1:2048,real(queso_calib_0_chain_inv_initPos0),'bo');
hold;
plot(1:2048,queso_calib_0_chain(1:2048,1),'r*');
plot(1:2048,real(queso_calib_0_chain_inv_initPos0)-queso_calib_0_chain(1:2048,1)','go'); %'
title('Difference between real component of uq_inv_fft and orginal information');

waitforbuttonpress;
clf;
plot(1:2048,imag(queso_calib_0_chain_inv_initPos0),'go');
title('Imag component of uq_inv_fft');
