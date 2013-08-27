sip1_filtChain_mh
sip2_filtChain_mh
sip3_filtChain_mh
sip4_filtChain_mh
sip5_filtChain_mh
dataPoints

b = 0.045213;
sigmaOrig = b/2;
n = [1 10 100 500 1000];

yMean = zeros(5,1);
yMean(1) = mean(yObs(1:n(1)));
yMean(2) = mean(yObs(1:n(2)));
yMean(3) = mean(yObs(1:n(3)));
yMean(4) = mean(yObs(1:n(4)));
yMean(5) = mean(yObs(1:n(5)));

qMean = yMean/b;

qStd = zeros(5,1);
qStd(1) = sigmaOrig/(sqrt(n(1)))/b;
qStd(2) = sigmaOrig/(sqrt(n(2)))/b;
qStd(3) = sigmaOrig/(sqrt(n(3)))/b;
qStd(4) = sigmaOrig/(sqrt(n(4)))/b;
qStd(5) = sigmaOrig/(sqrt(n(5)))/b;

%%%%%%%%%%%%%%%%%%%%

x1 = -0.2:0.001:0.2;
q1 = pdf('norm',x1,qMean(1),qStd(1));

[f,xi] = ksdensity(sip1_ip_mh_filtChain_unified(:,1),'function','pdf');
plot(x1,q1,'-r','linewidth',2)
hold
plot(xi,f,'--b','linewidth',2)
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1');
ylabel('Marginal posterior KDE');
title('2013-08-27 - Posterior - Verif. Prob. 5 - Case 1');
legend('Analytical',...
       'QUESO DRAM',...
       'location','northeast');
%a=axis();
%axis([-1 1 a(3) a(4)]);

print -dpng 2013_08_27_post_param01_verif_prob_5_case1.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

x2 = -0.1:0.001:0.1;
q2 = pdf('norm',x2,qMean(2),qStd(2));

[f,xi] = ksdensity(sip2_ip_mh_filtChain_unified(:,1),'function','pdf');
plot(x2,q2,'-r','linewidth',2)
hold
plot(xi,f,'--b','linewidth',2)
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1');
ylabel('Marginal posterior KDE');
title('2013-08-27 - Posterior - Verif. Prob. 5 - Case 2');
legend('Analytical',...
       'QUESO DRAM',...
       'location','northeast');
%a=axis();
%axis([-1 1 a(3) a(4)]);

print -dpng 2013_08_27_post_param01_verif_prob_5_case2.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

x3 = -0.02:0.0001:0.02;
q3 = pdf('norm',x3,qMean(3),qStd(3));

[f,xi] = ksdensity(sip3_ip_mh_filtChain_unified(:,1),'function','pdf');
plot(x3,q3,'-r','linewidth',2)
hold
plot(xi,f,'--b','linewidth',2)
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1');
ylabel('Marginal posterior KDE');
title('2013-08-27 - Posterior - Verif. Prob. 5 - Case 3');
legend('Analytical',...
       'QUESO DRAM',...
       'location','northeast');
%a=axis();
%axis([-1 1 a(3) a(4)]);

print -dpng 2013_08_27_post_param01_verif_prob_5_case3.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

x4 = -0.01:0.0001:0.01;
q4 = pdf('norm',x4,qMean(4),qStd(4));

[f,xi] = ksdensity(sip4_ip_mh_filtChain_unified(:,1),'function','pdf');
plot(x4,q4,'-r','linewidth',2)
hold
plot(xi,f,'--b','linewidth',2)
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1');
ylabel('Marginal posterior KDE');
title('2013-08-27 - Posterior - Verif. Prob. 5 - Case 4');
legend('Analytical',...
       'QUESO DRAM',...
       'location','northeast');
%a=axis();
%axis([-1 1 a(3) a(4)]);

print -dpng 2013_08_27_post_param01_verif_prob_5_case4.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

x5 = -0.002:0.00001:0.002;
q5 = pdf('norm',x5,qMean(5),qStd(5));

[f,xi] = ksdensity(sip5_ip_mh_filtChain_unified(:,1),'function','pdf');
plot(x5,q5,'-r','linewidth',2)
hold
plot(xi,f,'--b','linewidth',2)
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1');
ylabel('Marginal posterior KDE');
title('2013-08-27 - Posterior - Verif. Prob. 5 - Case 5');
legend('Analytical',...
       'QUESO DRAM',...
       'location','northeast');
%a=axis();
%axis([-1 1 a(3) a(4)]);

print -dpng 2013_08_27_post_param01_verif_prob_5_case5.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%
