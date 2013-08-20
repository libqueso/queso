cd ../tmp
load('pout.mat');
chain = zeros(10000,1);

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamWOs;
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1 (lambda\_eta)');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param01.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamUz(1);
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 2 (lambda\_w\_1)');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param02.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamUz(2);
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 3 (lambda\_w\_2)');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param03.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaU(1));
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. 1. avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 4 (rho\_w\_{1,1})');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param04.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaU(2));
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. 1. avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 5 (rho\_w\_{1,2})');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param05.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaU(3));
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. 1. avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 6 (rho\_w\_{2,1})');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param06.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaU(4));
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. 1. avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 7 (rho\_w\_{2,2})');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param07.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamWs(1);
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 8 (lambda\_s\_1)');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param08.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamWs(2);
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 9 (lambda\_s\_2)');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param09.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamOs;
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 10 (lambda\_y)');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param10.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamVz;
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 11 (lambda\_v\_1)');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param11.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaV);
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. 1. avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 12 (rho\_v\_{1,1})');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param12.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).theta;
end
[f,xi] = ksdensity(chain(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%avals=axis;
%axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter theta (13)');
ylabel('Marginal posterior KDE (GPMSA)');
print -dpng gpmsa_post_13_param13.png
waitforbuttonpress;
clf;

