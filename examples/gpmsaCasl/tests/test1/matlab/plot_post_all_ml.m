printGpmsaMatlab = 0;
printAxisMinMax  = 1;

if printGpmsaMatlab
load('gpmsa_pout.mat');
chain = zeros(10000,1);
end

priorSeq
filtChain_ml

%%%%%%%%%%%%%%%%%%%%

meanAll_ml = mean(gcm_ml_last_filtChain_unified,1)

stdAll_ml = std(gcm_ml_last_filtChain_unified,0,1)

waitforbuttonpress;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamWOs;
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,1)) min(chain(:,1)) min(gcm_priorSeq_unified(:,1))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,1)) max(chain(:,1)) max(gcm_priorSeq_unified(:,1))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,1)) min(gcm_priorSeq_unified(:,1))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,1)) max(gcm_priorSeq_unified(:,1))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,1),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,1),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([100       150       avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 1 (lambda\_eta)');
ylabel('Marginal posterior KDE');
print -dpng post_19_param01_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamUz(1);
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,2)) min(chain(:,1)) min(gcm_priorSeq_unified(:,2))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,2)) max(chain(:,1)) max(gcm_priorSeq_unified(:,2))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,2)) min(gcm_priorSeq_unified(:,2))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,2)) max(gcm_priorSeq_unified(:,2))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,2),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,2),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.9       1.1       avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 2 (lambda\_w\_1)');
ylabel('Marginal posterior KDE');
print -dpng post_19_param02_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamUz(2);
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,3)) min(chain(:,1)) min(gcm_priorSeq_unified(:,3))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,3)) max(chain(:,1)) max(gcm_priorSeq_unified(:,3))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,3)) min(gcm_priorSeq_unified(:,3))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,3)) max(gcm_priorSeq_unified(:,3))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,3),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,3),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.9       1.1       avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 3 (lambda\_w\_2)');
ylabel('Marginal posterior KDE');
print -dpng post_19_param03_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaU(1));
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,4)) min(chain(:,1)) min(gcm_priorSeq_unified(:,4))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,4)) max(chain(:,1)) max(gcm_priorSeq_unified(:,4))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,4)) min(gcm_priorSeq_unified(:,4))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,4)) max(gcm_priorSeq_unified(:,4))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,4),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,4),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.99      1.01      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 4 (rho\_w\_{1,1})');
ylabel('Marginal posterior KDE');
print -dpng post_19_param04_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaU(2));
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,5)) min(chain(:,1)) min(gcm_priorSeq_unified(:,5))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,5)) max(chain(:,1)) max(gcm_priorSeq_unified(:,5))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,5)) min(gcm_priorSeq_unified(:,5))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,5)) max(gcm_priorSeq_unified(:,5))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,5),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,5),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.99      1.01      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 5 (rho\_w\_{1,2})');
ylabel('Marginal posterior KDE');
print -dpng post_19_param05_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaU(3));
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,6)) min(chain(:,1)) min(gcm_priorSeq_unified(:,6))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,6)) max(chain(:,1)) max(gcm_priorSeq_unified(:,6))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,6)) min(gcm_priorSeq_unified(:,6))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,6)) max(gcm_priorSeq_unified(:,6))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,6),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,6),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.99      1.01      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 6 (rho\_w\_{1,3})');
ylabel('Marginal posterior KDE');
print -dpng post_19_param06_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaU(4));
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,7)) min(chain(:,1)) min(gcm_priorSeq_unified(:,7))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,7)) max(chain(:,1)) max(gcm_priorSeq_unified(:,7))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,7)) min(gcm_priorSeq_unified(:,7))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,7)) max(gcm_priorSeq_unified(:,7))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,7),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,7),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.99      1.01      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 7 (rho\_w\_{1,4})');
ylabel('Marginal posterior KDE');
print -dpng post_19_param07_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamWs(1);
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,8)) min(chain(:,1)) min(gcm_priorSeq_unified(:,8))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,8)) max(chain(:,1)) max(gcm_priorSeq_unified(:,8))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,8)) min(gcm_priorSeq_unified(:,8))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,8)) max(gcm_priorSeq_unified(:,8))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,8),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,8),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.99      1.01      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 8 (rho\_w\_{2,1})');
ylabel('Marginal posterior KDE');
print -dpng post_19_param08_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamWs(2);
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,9)) min(chain(:,1)) min(gcm_priorSeq_unified(:,9))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,9)) max(chain(:,1)) max(gcm_priorSeq_unified(:,9))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,9)) min(gcm_priorSeq_unified(:,9))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,9)) max(gcm_priorSeq_unified(:,9))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,9),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,9),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.99      1.01      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 9 (rho\_w\_{2,2})');
ylabel('Marginal posterior KDE');
print -dpng post_19_param09_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamOs;
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,10)) min(chain(:,1)) min(gcm_priorSeq_unified(:,10))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,10)) max(chain(:,1)) max(gcm_priorSeq_unified(:,10))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,10)) min(gcm_priorSeq_unified(:,10))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,10)) max(gcm_priorSeq_unified(:,10))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,10),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,10),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.99      1.01      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 10 (rho\_w\_{2,3})');
ylabel('Marginal posterior KDE');
print -dpng post_19_param10_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamVz;
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,11)) min(chain(:,1)) min(gcm_priorSeq_unified(:,11))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,11)) max(chain(:,1)) max(gcm_priorSeq_unified(:,11))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,11)) min(gcm_priorSeq_unified(:,11))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,11)) max(gcm_priorSeq_unified(:,11))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,11),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,11),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.99      1.01      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 11 (rho\_w\_{2,4})');
ylabel('Marginal posterior KDE');
print -dpng post_19_param11_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaV);
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,12)) min(chain(:,1)) min(gcm_priorSeq_unified(:,12))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,12)) max(chain(:,1)) max(gcm_priorSeq_unified(:,12))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,12)) min(gcm_priorSeq_unified(:,12))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,12)) max(gcm_priorSeq_unified(:,12))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,12),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,12),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([500       1000      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 12 (\lambda\_s\_1)');
ylabel('Marginal posterior KDE');
print -dpng post_19_param12_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).theta;
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,13)) min(chain(:,1)) min(gcm_priorSeq_unified(:,13))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,13)) max(chain(:,1)) max(gcm_priorSeq_unified(:,13))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,13)) min(gcm_priorSeq_unified(:,13))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,13)) max(gcm_priorSeq_unified(:,13))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,13),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,13),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([both_xmin 500       avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 13 (\lambda\_s\_2)');
ylabel('Marginal posterior KDE');
print -dpng post_19_param13_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).theta;
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,14)) min(chain(:,1)) min(gcm_priorSeq_unified(:,14))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,14)) max(chain(:,1)) max(gcm_priorSeq_unified(:,14))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,14)) min(gcm_priorSeq_unified(:,14))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,14)) max(gcm_priorSeq_unified(:,14))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,14),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,14),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([1.02      1.03      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 14 (\lambda\_y)');
ylabel('Marginal posterior KDE');
print -dpng post_19_param14_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).theta;
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,15)) min(chain(:,1)) min(gcm_priorSeq_unified(:,15))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,15)) max(chain(:,1)) max(gcm_priorSeq_unified(:,15))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,15)) min(gcm_priorSeq_unified(:,15))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,15)) max(gcm_priorSeq_unified(:,15))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,15),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,15),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([1600      1700      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 15 (\lambda\_v\_1)');
ylabel('Marginal posterior KDE');
print -dpng post_19_param15_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).theta;
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,16)) min(chain(:,1)) min(gcm_priorSeq_unified(:,16))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,16)) max(chain(:,1)) max(gcm_priorSeq_unified(:,16))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,16)) min(gcm_priorSeq_unified(:,16))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,16)) max(gcm_priorSeq_unified(:,16))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,16),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,16),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.99      both_xmax avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 16 (\rho\_v\_1)');
ylabel('Marginal posterior KDE');
print -dpng post_19_param16_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).theta;
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,17)) min(chain(:,1)) min(gcm_priorSeq_unified(:,17))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,17)) max(chain(:,1)) max(gcm_priorSeq_unified(:,17))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,17)) min(gcm_priorSeq_unified(:,17))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,17)) max(gcm_priorSeq_unified(:,17))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,17),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,17),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.14      0.16      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 17 (\theta\_1)');
ylabel('Marginal posterior KDE');
print -dpng post_19_param17_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).theta;
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,18)) min(chain(:,1)) min(gcm_priorSeq_unified(:,18))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,18)) max(chain(:,1)) max(gcm_priorSeq_unified(:,18))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,18)) min(gcm_priorSeq_unified(:,18))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,18)) max(gcm_priorSeq_unified(:,18))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,18),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,18),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.42      0.44      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 18 (\theta\_2)');
ylabel('Marginal posterior KDE');
print -dpng post_19_param18_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

if printGpmsaMatlab
for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).theta;
end
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,19)) min(chain(:,1)) min(gcm_priorSeq_unified(:,19))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,19)) max(chain(:,1)) max(gcm_priorSeq_unified(:,19))]);
else
both_xmin = min([min(gcm_ml_last_filtChain_unified(:,19)) min(gcm_priorSeq_unified(:,19))]);
both_xmax = max([max(gcm_ml_last_filtChain_unified(:,19)) max(gcm_priorSeq_unified(:,19))]);
end

[h,hi] = ksdensity(gcm_priorSeq_unified(:,19),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold

if printGpmsaMatlab
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
end

[ffilt,xi] = ksdensity(gcm_ml_last_filtChain_unified(:,19),'function','pdf');
plot(xi,ffilt,'-b','linewidth',4);

avals=axis;
if printAxisMinMax
axis([both_xmin both_xmax avals(3) min([avals(4) avals(4)])]);
else
axis([0.70      0.72      avals(3) avals(4)]);
end
grid minor;
set(gca,'FontSize',16);
if printGpmsaMatlab
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO ML',...
       'location','north');
else
legend('Prior KDE',...
       'Post QUESO ML',...
       'location','north');
end
xlabel('Parameter 19 (\theta\_3)');
ylabel('Marginal posterior KDE');
print -dpng post_19_param19_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%
