cd outputData
sipOutput_sub0
sipExtraOutput_sub0
cd ..
load nada

d=nada(:,2);    
[nsimu,temp]=size(d);
[temp,npar]=size(ip_mc_chain_sub0);
meanVec=[0 0 0 0];
covMatrix = [0.263584 0.0820654 0.0748222 -0.410472 
0.0820654 0.0414596 0.0241198 -0.137645 
0.0748222 0.0241198 0.0511172 -0.140059 
-0.410472 -0.137645 -0.140059 0.698176 
];

path(path,'/h2/prudenci/DownloadedPackages/mcmc02Jun2008/code');    
c50  = chiqf(0.50,npar);
c95  = chiqf(0.95,npar);
cc50 = sum(d<c50)./nsimu;
cc95 = sum(d<c95)./nsimu;

plot(ip_mc_chain_sub0(:,1),ip_mc_chain_sub0(:,2),'.');
xlabel(ip_mc_componentNames(1));
ylabel(ip_mc_componentNames(2));
title(sprintf('Rejected = %.1f%%, c50 = %.1f%%, c95 = %.1f%%', ...
              ip_mc_rejected*100, cc50*100, cc95*100))
hold

c50  = chiqf(0.50,2);
c95  = chiqf(0.95,2);
ellipse(meanVec(1:2),c50*covMatrix(1:2,1:2),'r--','LineWidth',2);
ellipse(meanVec(1:2),c95*covMatrix(1:2,1:2),'r--','LineWidth',2);
waitforbuttonpress;
clf;

plot(ip_mc_chain_GkdePosits_sub0(1,:),ip_mc_chain_GkdeValues_sub0(1,:),'-b');
hold
plot(ip_mc_chain_GkdePosits_sub0(2,:),ip_mc_chain_GkdeValues_sub0(2,:),'-r');
plot(ip_mc_chain_GkdePosits_sub0(3,:),ip_mc_chain_GkdeValues_sub0(3,:),'--b');
plot(ip_mc_chain_GkdePosits_sub0(4,:),ip_mc_chain_GkdeValues_sub0(4,:),'--r');

ylabel('marginal pdf','fontsize',20);
xlabel('Parameter values','fontsize',20);
title('Posterior marginal pdfs of parameters','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('param1',...
       'param2',...
       'param3',...
       'param4',...
       'location','southwest');
