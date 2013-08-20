path(path,'/Users/prudenci/works/personal/svns/sa/stalg/gp_lanl/code_orig/gpmsa_Sourceforge_0609_changed_locally/matlab');
load('gpmsa_pout.mat');
%cd ../outputData_m1_2012_12_29__02_08_hs_CT__lonestar
%cd ../outputData_m1_2012_12_30__17_40_hs_CT__mac
cd ../outputData_m1_2013_05_12__03_19_hs_CT__mac
predictionGrid1_sub0
predictionGrid2_sub0
w1mat_sub0
w2mat_sub0
cd ../matlab

%%%%%%%%%%%%%%%%%%%%

surf(towerPredictionGrid1_sub0,towerPredictionGrid2_sub0,towerW1mat_sub0'); %'
a = axis;
axis([0 1 0 1 a(5) a(6)]);
grid minor;
set(gca,'FontSize',16);
xlabel('x');
ylabel('theta');
zlabel('Mean(w_1) QUESO');
print -dpng pred_w_1_queso.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%

surf(towerPredictionGrid1_sub0,towerPredictionGrid2_sub0,towerW2mat_sub0'); %'
a = axis;
axis([0 1 0 1 a(5) a(6)]);
grid minor;
set(gca,'FontSize',16);
xlabel('x');
ylabel('theta');
zlabel('Mean(w_2) QUESO');
print -dpng pred_w_2_queso.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%

from = 2000; % ernesto 200 2000
to = length(pout.pvals);
thismany = 500; % ernesto 100 500
ip = round(linspace(from,to,thismany));

model = pout.model;
data = pout.data;
pvals = pout.pvals(ip);
gridNumPos = 16;
gridDelta = 1./(gridNumPos-1);
gridUni = 0:gridDelta:1;

[gridx, gridy] = meshgrid(gridUni, gridUni);
npc = size(pout.simData.Ksim, 2); % number of principal components
xpred = gridx(:); 
theta = gridy(:);
pred = gPred(xpred, pvals, model, data, 'wpred', theta);
pm = squeeze(mean(pred.w, 1));

%%%%%%%%%%%%%%%%%%%%

w1 = reshape(pm(1, :), [gridNumPos gridNumPos]); % make the mean for component ii match the grid size
surf(gridUni',gridUni',w1);
a = axis;
axis([0 1 0 1 a(5) a(6)]);
grid minor;
set(gca,'FontSize',16);
xlabel('x');
ylabel('theta');
zlabel('Mean(w_1) GPMSA');
title('PC 1');
print -dpng pred_w_1_gpmsa.png

waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

w2 = reshape(pm(2, :), [gridNumPos gridNumPos]); % make the mean for component ii match the grid size
surf(gridUni',gridUni',w2);
a = axis;
axis([0 1 0 1 a(5) a(6)]);
grid minor;
set(gca,'FontSize',16);
xlabel('x');
ylabel('\theta');
zlabel('Mean(w_2) GPMSA');
title('PC 2');
print -dpng pred_w_2_gpmsa.png

waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

surf(towerPredictionGrid1_sub0,towerPredictionGrid2_sub0,towerW1mat_sub0'-w1); %'
a = axis;
axis([0 1 0 1 a(5) a(6)]);
grid minor;
set(gca,'FontSize',16);
xlabel('x');
ylabel('theta');
zlabel('Mean(w_1) Diff');
print -dpng pred_w_1_diff.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%

surf(towerPredictionGrid1_sub0,towerPredictionGrid2_sub0,towerW2mat_sub0'-w2); %'
a = axis;
axis([0 1 0 1 a(5) a(6)]);
grid minor;
set(gca,'FontSize',16);
xlabel('x');
ylabel('theta');
zlabel('Mean(w_2) Diff');
print -dpng pred_w_2_diff.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%

surf(towerPredictionGrid1_sub0,towerPredictionGrid2_sub0,log10(abs((w1-towerW1mat_sub0')./w1))); %'
a = axis;
axis([0 1 0 1 a(5) a(6)]);
grid minor;
set(gca,'FontSize',16);
xlabel('x');
ylabel('theta');
zlabel('Mean(w_1) log10(|Diff/GPMSA|)');
print -dpng pred_w_1_diff_rel.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%

surf(towerPredictionGrid1_sub0,towerPredictionGrid2_sub0,log10(abs((w2-towerW2mat_sub0')./w2))); %'
a = axis;
axis([0 1 0 1 a(5) a(6)]);
grid minor;
set(gca,'FontSize',16);
xlabel('x');
ylabel('theta');
zlabel('Mean(w_2) log10(|Diff/GPMSA|)');
print -dpng pred_w_2_diff_rel.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%

