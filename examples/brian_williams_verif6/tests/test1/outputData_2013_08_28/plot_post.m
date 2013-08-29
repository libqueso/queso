filtChain_mh
dataPoints

%%%%%%%%%%%%%%%%%%%%

dat = textread('ctf-2p.txt');
resid = dat(:,2)-17333.65*dat(:,1)+9378.29*dat(:,1).^2-121370.92;
ISigma = diag(kron([4.148751 4.170606 4.190899 4.209629],ones(1,15)).^(-2));
X = ones(60,2);
X(:,1) = 50.45;
X(:,2) = -2.99;
PCov = (diag([.1/3 .1/3].^(-2))+X'*ISigma*X)\eye(2); %'
PMean = PCov*(X'*ISigma*resid); %'

PMean = PMean

PCov = PCov

%%%%%%%%%%%%%%%%%%%%

x1 = -0.06:0.0001:0.06;
q1 = pdf('norm',x1,PMean(1),sqrt(PCov(1,1)));
[f,xi] = ksdensity(sip_ip_mh_filtChain_unified(:,1),'function','pdf');
plot(x1,q1,'-r','linewidth',2)
hold
plot(xi,f,'--b','linewidth',2)
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1');
ylabel('Marginal posterior KDE');
title('2013-08-28 - Posterior - Verif. Prob. 6');
legend('Analytical',...
       'QUESO DRAM',...
       'location','northeast');
%a=axis();
%axis([-1 1 a(3) a(4)]);

print -dpng 2013_08_28_post_param01_verif_prob_6.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

x2 = -0.2:0.001:0.2;
q2 = pdf('norm',x2,PMean(2),sqrt(PCov(2,2)));
[f,xi] = ksdensity(sip_ip_mh_filtChain_unified(:,2),'function','pdf');
plot(x2,q2,'-r','linewidth',2)
hold
plot(xi,f,'--b','linewidth',2)
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 2');
ylabel('Marginal posterior KDE');
title('2013-08-28 - Posterior - Verif. Prob. 6');
legend('Analytical',...
       'QUESO DRAM',...
       'location','northeast');
%a=axis();
%axis([-1 1 a(3) a(4)]);

print -dpng 2013_08_28_post_param02_verif_prob_6.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

x1 = -0.048:0.001:0.068;
x2 = -0.185:0.001:0.185;

[m1,n1]=size(x1);
[m2,n2]=size(x2);
yMat=zeros(n1,n2);

tmp = 1./2/pi/sqrt(det(PCov));
for i=1:n1
for j=1:n2
  xVec    = [x1(i) x2(j)]'; %'
  diffVec = xVec - PMean;
  tmpVec = PCov\diffVec;
  yMat(i,j) = tmp * exp(-0.5*diffVec'*tmpVec); %'
end
end

surf(x1,x2,yMat'); %'
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1');
ylabel('Parameter 2');
zlabel('Analytical 2D posterior KDE');
%a=axis();
%axis([5 40 4 22 a(5) a(6)]);

view([315,45]);
title('2013-08-28 - ANALYTICAL Posterior - Verif. Prob. 6 - View A');
print -dpng 2013_08_28_analytical_post_kde_params_1_and_2_verif_prob6_viewA.png
waitforbuttonpress;

view([225,45]);
title('2013-08-28 - ANALYTICAL Posterior - Verif. Prob. 6 - View B');
print -dpng 2013_08_28_analytical_post_kde_params_1_and_2_verif_prob6_viewB.png
waitforbuttonpress;

view([135,45]);
title('2013-08-28 - ANALYTICAL Posterior - Verif. Prob. 6 - View C');
print -dpng 2013_08_28_analytical_post_kde_params_1_and_2_verif_prob6_viewC.png
waitforbuttonpress;

view([45,45]);
title('2013-08-28 - ANALYTICAL Posterior - Verif. Prob. 6 - View D');
print -dpng 2013_08_28_analytical_post_kde_params_1_and_2_verif_prob6_viewD.png
waitforbuttonpress;

view([0,90]);
title('2013-08-28 - ANALYTICAL Posterior - Verif. Prob. 6 - View E');
print -dpng 2013_08_28_analytical_post_kde_params_1_and_2_verif_prob6_viewE.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

data=[sip_ip_mh_filtChain_unified(:,1), sip_ip_mh_filtChain_unified(:,2)];
[bandwidth,density,X,Y]=kde2d(data,64);
surf(X,Y,density);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1');
ylabel('Parameter 2');
zlabel('QUESO 2D posterior KDE');
%a=axis();
%axis([5 40 4 22 a(5) a(6)]);

view([315,45]);
title('2013-08-28 - QUESO Posterior - Verif. Prob. 6 - View A');
print -dpng 2013_08_28_queso_post_kde_params_1_and_2_verif_prob6_viewA.png
waitforbuttonpress;

view([225,45]);
title('2013-08-28 - QUESO Posterior - Verif. Prob. 6 - View B');
print -dpng 2013_08_28_queso_post_kde_params_1_and_2_verif_prob6_viewB.png
waitforbuttonpress;

view([135,45]);
title('2013-08-28 - QUESO Posterior - Verif. Prob. 6 - View C');
print -dpng 2013_08_28_queso_post_kde_params_1_and_2_verif_prob6_viewC.png
waitforbuttonpress;

view([45,45]);
title('2013-08-28 - QUESO Posterior - Verif. Prob. 6 - View D');
print -dpng 2013_08_28_queso_post_kde_params_1_and_2_verif_prob6_viewD.png
waitforbuttonpress;

view([0,90]);
title('2013-08-28 - QUESO Posterior - Verif. Prob. 6 - View E');
print -dpng 2013_08_28_queso_post_kde_params_1_and_2_verif_prob6_viewE.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%
