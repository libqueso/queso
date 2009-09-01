cd outputData
sipOutput_sub0
file_sip_raw
cd ..
load nada

d=nada(:,1);    
[nsimu,temp]=size(d);
[temp,npar]=size(ip_mc_rawChain_unified);
meanVec=[0 0 0 0];
covMatrix = [0.263584 0.0820654 0.0748222 -0.410472 
0.0820654 0.0414596 0.0241198 -0.137645 
0.0748222 0.0241198 0.0511172 -0.140059 
-0.410472 -0.137645 -0.140059 0.698176 
];

c50  = chi2inv(0.50,npar);
c95  = chi2inv(0.95,npar);
cc50 = sum(d<c50)./nsimu;
cc95 = sum(d<c95)./nsimu;

plot(ip_mc_rawChain_unified(:,1),ip_mc_rawChain_unified(:,2),'.');
xlabel('\theta_1');
ylabel('\theta_2');
title(sprintf('Rejected = %.1f%%, c50 = %.1f%%, c95 = %.1f%%', ...
              ip_mc_rejected*100, cc50*100, cc95*100))
hold

c50  = chi2inv(0.50,2);
c95  = chi2inv(0.95,2);

t = linspace(0,2*pi)'; %'
R = chol(c50*covMatrix(1:2,1:2));
x = meanVec(1) + R(1,1).*cos(t);
y = meanVec(2) + R(1,2).*cos(t) + R(2,2).*sin(t);
plot(x,y,'r--','LineWidth',2);

R = chol(c95*covMatrix(1:2,1:2));
x = meanVec(1) + R(1,1).*cos(t);
y = meanVec(2) + R(1,2).*cos(t) + R(2,2).*sin(t);
plot(x,y,'r--','LineWidth',2);

print -dpng fig1.png
waitforbuttonpress;
clf;

plot(ip_mc_rawChain_unifGkdePosits_sub0(1,:),ip_mc_rawChain_unifGkdeValues_sub0(1,:),'-b');
hold
plot(ip_mc_rawChain_unifGkdePosits_sub0(2,:),ip_mc_rawChain_unifGkdeValues_sub0(2,:),'-r');
plot(ip_mc_rawChain_unifGkdePosits_sub0(3,:),ip_mc_rawChain_unifGkdeValues_sub0(3,:),'--b');
plot(ip_mc_rawChain_unifGkdePosits_sub0(4,:),ip_mc_rawChain_unifGkdeValues_sub0(4,:),'--r');

ylabel('marginal pdf','fontsize',20);
xlabel('Parameter values','fontsize',20);
title('Posterior marginal pdfs of parameters','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('param1',...
       'param2',...
       'param3',...
       'param4',...
       'location','northwest');
print -dpng fig2.png
