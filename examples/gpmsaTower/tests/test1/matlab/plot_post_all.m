printRawChain = 0;

load('gpmsa_pout.mat');
chain = zeros(10000,1);
%cd ../outputData_m1_2012_12_29__02_08_hs_CT__lonestar
%cd ../outputData_m1_2012_12_30__17_40_hs_CT__mac
cd ../outputData_m1_2013_05_12__03_19_hs_CT__mac
priorSeq
rawChain_mh
filtChain_mh
cd ../matlab

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamWOs;
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,1),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,1),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,1),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,1)) min(gcm_mh_filtChain_unified(:,1)) min(chain(:,1)) min(gcm_priorSeq_unified(:,1))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,1)) max(gcm_mh_filtChain_unified(:,1)) max(chain(:,1)) max(gcm_priorSeq_unified(:,1))]);
avals=axis;
axis([both_xmin both_xmax avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 1 (lambda\_eta)');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param01_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamUz(1);
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,2),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,2),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,2),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,2)) min(gcm_mh_filtChain_unified(:,2)) min(chain(:,1)) min(gcm_priorSeq_unified(:,2))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,2)) max(gcm_mh_filtChain_unified(:,2)) max(chain(:,1)) max(gcm_priorSeq_unified(:,2))]);
avals=axis;
axis([both_xmin both_xmax avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 2 (lambda\_w\_1)');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param02_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamUz(2);
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,3),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,3),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,3),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,3)) min(gcm_mh_filtChain_unified(:,3)) min(chain(:,1)) min(gcm_priorSeq_unified(:,3))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,3)) max(gcm_mh_filtChain_unified(:,3)) max(chain(:,1)) max(gcm_priorSeq_unified(:,3))]);
avals=axis;
axis([both_xmin both_xmax avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 3 (lambda\_w\_2)');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param03_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaU(1));
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,4),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,4),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,4),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,4)) min(gcm_mh_filtChain_unified(:,4)) min(chain(:,1)) min(gcm_priorSeq_unified(:,4))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,4)) max(gcm_mh_filtChain_unified(:,4)) max(chain(:,1)) max(gcm_priorSeq_unified(:,4))]);
avals=axis;
axis([both_xmin both_xmax avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 4 (rho\_w\_{1,1})');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param04_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaU(2));
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,5),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,5),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,5),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,5)) min(gcm_mh_filtChain_unified(:,5)) min(chain(:,1)) min(gcm_priorSeq_unified(:,5))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,5)) max(gcm_mh_filtChain_unified(:,5)) max(chain(:,1)) max(gcm_priorSeq_unified(:,5))]);
avals=axis;
axis([both_xmin both_xmax avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 5 (rho\_w\_{1,2})');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param05_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaU(3));
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,6),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,6),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,6),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,6)) min(gcm_mh_filtChain_unified(:,6)) min(chain(:,1)) min(gcm_priorSeq_unified(:,6))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,6)) max(gcm_mh_filtChain_unified(:,6)) max(chain(:,1)) max(gcm_priorSeq_unified(:,6))]);
avals=axis;
axis([both_xmin both_xmax avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 6 (rho\_w\_{2,1})');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param06_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaU(4));
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,7),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,7),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,7),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,7)) min(gcm_mh_filtChain_unified(:,7)) min(chain(:,1)) min(gcm_priorSeq_unified(:,7))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,7)) max(gcm_mh_filtChain_unified(:,7)) max(chain(:,1)) max(gcm_priorSeq_unified(:,7))]);
avals=axis;
axis([both_xmin both_xmax avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 7 (rho\_w\_{2,2})');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param07_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamWs(1);
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,8),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,8),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,8),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,8)) min(gcm_mh_filtChain_unified(:,8)) min(chain(:,1)) min(gcm_priorSeq_unified(:,8))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,8)) max(gcm_mh_filtChain_unified(:,8)) max(chain(:,1)) max(gcm_priorSeq_unified(:,8))]);
avals=axis;
axis([both_xmin both_xmax avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 8 (lambda\_s\_1)');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param08_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamWs(2);
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,9),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,9),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,9),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,9)) min(gcm_mh_filtChain_unified(:,9)) min(chain(:,1)) min(gcm_priorSeq_unified(:,9))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,9)) max(gcm_mh_filtChain_unified(:,9)) max(chain(:,1)) max(gcm_priorSeq_unified(:,9))]);
avals=axis;
axis([both_xmin both_xmax avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 9 (lambda\_s\_2)');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param09_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamOs;
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,10),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,10),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,10),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,10)) min(gcm_mh_filtChain_unified(:,10)) min(chain(:,1)) min(gcm_priorSeq_unified(:,10))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,10)) max(gcm_mh_filtChain_unified(:,10)) max(chain(:,1)) max(gcm_priorSeq_unified(:,10))]);
avals=axis;
axis([both_xmin both_xmax avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 10 (lambda\_y)');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param10_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).lamVz;
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,11),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,11),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,11),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,11)) min(gcm_mh_filtChain_unified(:,11)) min(chain(:,1)) min(gcm_priorSeq_unified(:,11))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,11)) max(gcm_mh_filtChain_unified(:,11)) max(chain(:,1)) max(gcm_priorSeq_unified(:,11))]);
avals=axis;
axis([both_xmin 2000 avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 11 (lambda\_v\_1)');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param11_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = exp(-0.25*pout.pvals(1,1300+i).betaV);
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,12),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,12),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,12),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,12)) min(gcm_mh_filtChain_unified(:,12)) min(chain(:,1)) min(gcm_priorSeq_unified(:,12))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,12)) max(gcm_mh_filtChain_unified(:,12)) max(chain(:,1)) max(gcm_priorSeq_unified(:,12))]);
avals=axis;
axis([both_xmin both_xmax avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 12 (rho\_v\_{1,1})');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param12_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

for i=1:10000
  chain(i,1) = pout.pvals(1,1300+i).theta;
end
[h,hi] = ksdensity(gcm_priorSeq_unified(:,13),'function','pdf');
plot(hi,h,'-g','linewidth',2)
hold
[g,gi] = ksdensity(chain(:,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
if printRawChain
[fraw,xi] = ksdensity(gcm_mh_rawChain_unified(:,13),'function','pdf');
plot(xi,fraw,':b','linewidth',2);
end
[ffilt,xi] = ksdensity(gcm_mh_filtChain_unified(:,13),'function','pdf');
plot(xi,ffilt,'-b','linewidth',2);
both_xmin = min([min(gcm_mh_rawChain_unified(:,13)) min(gcm_mh_filtChain_unified(:,13)) min(chain(:,1)) min(gcm_priorSeq_unified(:,13))]);
both_xmax = max([max(gcm_mh_rawChain_unified(:,13)) max(gcm_mh_filtChain_unified(:,13)) max(chain(:,1)) max(gcm_priorSeq_unified(:,13))]);
avals=axis;
axis([both_xmin both_xmax avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
if printRawChain
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO raw',...
       'Post QUESO filt',...
       'location','north');
else
legend('Prior KDE',...
       'Post GPMSA',...
       'Post QUESO',...
       'location','north');
end
xlabel('Parameter 13 (theta)');
ylabel('Marginal posterior KDE (BOTH + filter)');
print -dpng post_13_param13_all.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%
