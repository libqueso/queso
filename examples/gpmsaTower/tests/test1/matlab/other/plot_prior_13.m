cd ../outputData_m1_2012_12_12__11_11_hs_CT
priorSeq

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,1),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,1));
gcm_xmax = max(gcm_priorSeq_unified(:,1));
gcm_xdelta = (gcm_xmax - gcm_xmin)/100;
gcm_ymax = max(h);
a=5;
b=0.005;
x=gcm_xmin:gcm_xdelta:gcm_xmax;
normalizationFactor = b^a/gamma(a);
y=normalizationFactor*x.^(a-1).*exp(-b.*x);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
avals=axis;
axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 1 (lambda\_eta)');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param01.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,2),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,2));
gcm_xmax = max(gcm_priorSeq_unified(:,2));
gcm_xdelta = (gcm_xmax - gcm_xmin)/100;
gcm_ymax = max(h);
a=5;
b=5;
x=gcm_xmin:gcm_xdelta:gcm_xmax;
normalizationFactor = b^a/gamma(a);
y=normalizationFactor*x.^(a-1).*exp(-b.*x);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
avals=axis;
axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 2 (lambda\_w\_1)');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param02.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,3),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,3));
gcm_xmax = max(gcm_priorSeq_unified(:,3));
gcm_xdelta = (gcm_xmax - gcm_xmin)/100;
gcm_ymax = max(h);
a=5;
b=5;
x=gcm_xmin:gcm_xdelta:gcm_xmax;
normalizationFactor = b^a/gamma(a);
y=normalizationFactor*x.^(a-1).*exp(-b.*x);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
avals=axis;
axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 3 (lambda\_w\_2)');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param03.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,4),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,4));
gcm_xmax = max(gcm_priorSeq_unified(:,4));
gcm_xdelta = (gcm_xmax - gcm_xmin)/100;
gcm_ymax = max(h);
a=1;
b=0.1;
x=gcm_xmin:gcm_xdelta:0.99;
normalizationFactor = 1./beta(a,b);
y=normalizationFactor*x.^(a-1).*(1-x).^(b-1);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
avals=axis;
axis([0. 1. avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 4 (rho\_w\_{1,1})');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param04.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,5),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,5));
gcm_xmax = max(gcm_priorSeq_unified(:,5));
gcm_xdelta = (gcm_xmax - gcm_xmin)/1000;
gcm_ymax = max(h);
a=1;
b=0.1;
x=gcm_xmin:gcm_xdelta:0.99;
normalizationFactor = 1./beta(a,b);
y=normalizationFactor*x.^(a-1).*(1-x).^(b-1);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
avals=axis;
axis([0. 1. avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 5 (rho\_w\_{1,2})');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param05.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,6),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,6));
gcm_xmax = max(gcm_priorSeq_unified(:,6));
gcm_xdelta = (gcm_xmax - gcm_xmin)/1000;
gcm_ymax = max(h);
a=1;
b=0.1;
x=gcm_xmin:gcm_xdelta:0.99;
normalizationFactor = 1./beta(a,b);
y=normalizationFactor*x.^(a-1).*(1-x).^(b-1);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
avals=axis;
axis([0. 1. avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 6 (rho\_w\_{2,1})');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param06.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,7),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,7));
gcm_xmax = max(gcm_priorSeq_unified(:,7));
gcm_xdelta = (gcm_xmax - gcm_xmin)/1000;
gcm_ymax = max(h);
a=1;
b=0.1;
x=gcm_xmin:gcm_xdelta:0.99;
normalizationFactor = 1./beta(a,b);
y=normalizationFactor*x.^(a-1).*(1-x).^(b-1);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
avals=axis;
axis([0. 1. avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 7 (rho\_w\_{2,2})');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param07.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,8),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,8));
gcm_xmax = max(gcm_priorSeq_unified(:,8));
gcm_xdelta = (gcm_xmax - gcm_xmin)/100;
gcm_ymax = max(h);
a=3;
b=0.003;
x=gcm_xmin:gcm_xdelta:gcm_xmax;
normalizationFactor = b^a/gamma(a);
y=normalizationFactor*x.^(a-1).*exp(-b.*x);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
avals=axis;
axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 8 (lambda\_s\_1)');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param08.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,9),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,9));
gcm_xmax = max(gcm_priorSeq_unified(:,9));
gcm_xdelta = (gcm_xmax - gcm_xmin)/100;
gcm_ymax = max(h);
a=3;
b=0.003;
x=gcm_xmin:gcm_xdelta:gcm_xmax;
normalizationFactor = b^a/gamma(a);
y=normalizationFactor*x.^(a-1).*exp(-b.*x);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
avals=axis;
axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 9 (lambda\_s\_2)');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param09.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,10),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,10));
gcm_xmax = max(gcm_priorSeq_unified(:,10));
gcm_xdelta = (gcm_xmax - gcm_xmin)/100;
gcm_ymax = max(h);
a=1;
b=0.001;
x=gcm_xmin:gcm_xdelta:gcm_xmax;
normalizationFactor = b^a/gamma(a);
y=normalizationFactor*x.^(a-1).*exp(-b.*x);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
avals=axis;
axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 10 (lambda\_y)');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param10.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,11),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,11));
gcm_xmax = max(gcm_priorSeq_unified(:,11));
gcm_xdelta = (gcm_xmax - gcm_xmin)/100;
gcm_ymax = max(h);
a=1;
b=0.001;
x=gcm_xmin:gcm_xdelta:gcm_xmax;
normalizationFactor = b^a/gamma(a);
y=normalizationFactor*x.^(a-1).*exp(-b.*x);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
avals=axis;
axis([0. avals(2) avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 11 (lambda\_v\_1)');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param11.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,12),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,12));
gcm_xmax = max(gcm_priorSeq_unified(:,12));
gcm_xdelta = (gcm_xmax - gcm_xmin)/1000;
gcm_ymax = max(h);
a=1;
b=0.1;
x=gcm_xmin:gcm_xdelta:0.99;
normalizationFactor = 1./beta(a,b);
y=normalizationFactor*x.^(a-1).*(1-x).^(b-1);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
avals=axis;
axis([0. 1. avals(3) avals(4)]);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 12 (rho\_v\_{1,1})');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param12.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[h,hi] = ksdensity(gcm_priorSeq_unified(:,13),'function','pdf');
plot(hi,h,'-g','linewidth',2)
%waitforbuttonpress;
hold
gcm_xmin = min(gcm_priorSeq_unified(:,13));
gcm_xmax = max(gcm_priorSeq_unified(:,13));
gcm_xdelta = (gcm_xmax - gcm_xmin)/1000;
gcm_ymax = max(h);
m=0.5;
v=100.;
x=gcm_xmin:gcm_xdelta:gcm_xmax;
normalizationFactor=1.; %/sqrt(2*pi*v);
y=normalizationFactor*exp(-0.5*(x-m).*(x-m)/v);
yfactor=1; %gcm_ymax/max(y);
y=y.*yfactor;
plot(x,y,'--k','linewidth',2);
grid minor;
set(gca,'FontSize',16);
legend('Prior KDE',...
       'Exact Pdf',...
       'location','north');
xlabel('Parameter 13 (theta)');
ylabel('Marginal prior KDE');
print -dpng queso_prior_13_param13.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%
