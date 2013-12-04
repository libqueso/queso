%cd ../test/test_2013_12_02/outputData
cd outputData

x=0.:.1:150.;

v1=1;
v2=25;

y1=(x-10).*(x-10)/(2*v1);
y2=(x-100).*(x-100)/(2*v2);

f1=(exp(-y1))/(sqrt(v1)*sqrt(2*pi));
f2=(exp(-y2))/(sqrt(v2)*sqrt(2*pi));

f=f1/2+f2/2;
plot(x,f,'k-','linewidth',2);
ylabel('f(\theta)','fontsize',16);
xlabel('\theta','fontsize',16);
grid minor
set(gca,'fontsize',16);
title('Bimodal likelihood function','fontsize',20);

print -dpng bimodal_likelihood.png

%cd ../../../matlab
cd ..
