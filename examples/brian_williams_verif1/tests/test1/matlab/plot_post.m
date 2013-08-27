filtChain_mh
dataPoints

a = 126831.7;
b = 112136.1;
sigmaOrig = 4229.55;
[n,m] = size(yObs);
yMean = mean(yObs);
qMean = (yMean - a)/b;
qStd = sigmaOrig/(sqrt(n))/b;

x = -0.01:0.0001:0.01;
q = pdf('norm',x,qMean,qStd);

%%%%%%%%%%%%%%%%%%%%

[f,xi] = ksdensity(sip_ip_mh_filtChain_unified(:,1),'function','pdf');
plot(x,q,'-r','linewidth',2)
hold
plot(xi,f,'--b','linewidth',2)
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1');
ylabel('Marginal posterior KDE');
title('2013-08-26 - Posterior - Verif. Prob. 1');
legend('Analytical',...
       'QUESO DRAM',...
       'location','northeast');
%a=axis();
%axis([-1 1 a(3) a(4)]);

print -dpng 2013_08_26_post_param01_verif_prob_1.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%
