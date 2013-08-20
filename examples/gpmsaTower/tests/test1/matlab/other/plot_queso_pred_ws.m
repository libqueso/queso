cd ../outputData_m1_2012_12_25__21_38_hs_CT__mac
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
zlabel('Mean(w_1)');
print -dpng queso_pred_w_1.png
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
zlabel('Mean(w_2)');
print -dpng queso_pred_w_2.png
waitforbuttonpress;
clf
