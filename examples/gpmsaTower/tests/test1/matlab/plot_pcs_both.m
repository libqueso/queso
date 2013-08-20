load('gpmsa_w.txt');
%cd ../outputData_m1_2012_12_29__02_08_hs_CT__lonestar
%cd ../outputData_m1_2012_12_30__17_40_hs_CT__mac
cd ../outputData_m1_2013_05_12__03_19_hs_CT__mac
mat_K_eta_sub0
cd ../matlab

[m n] = size(Kmat_eta_sub0);

%%%%%%%%%%%%%%%%%%%%

plot(1:m,gpmsa_w(:,1),'-ro','linewidth',2);
hold
plot(1:m,Kmat_eta_sub0(:,1),':bx','linewidth',2);
grid minor;
a = axis;
axis([1 m a(3) a(4)]);
set(gca,'FontSize',16);
legend('GPMSA',...
       'QUESO',...
       'location','north');
xlabel('Index');
ylabel('Simulation PC 1');
print -dpng simul_pc_1_both.png
waitforbuttonpress;
clf

plot(1:m,gpmsa_w(:,1)-Kmat_eta_sub0(:,1),'-bo','linewidth',2);
grid minor;
a = axis;
axis([1 m a(3) a(4)]);
set(gca,'FontSize',16);
xlabel('Index');
ylabel('Sim. PC 1, Diff');
print -dpng simul_pc_1_diff.png
waitforbuttonpress;
clf

plot(1:m,log10(abs((gpmsa_w(:,1)-Kmat_eta_sub0(:,1))./gpmsa_w(:,1))),'-bo','linewidth',2);
grid minor;
a = axis;
axis([1 m a(3) a(4)]);
set(gca,'FontSize',16);
xlabel('Index');
ylabel('Sim. PC 1, log10(|Diff/GPMSA|)');
print -dpng simul_pc_1_diff_rel.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%

plot(1:m,gpmsa_w(:,2),'-ro','linewidth',2);
hold
plot(1:m,Kmat_eta_sub0(:,2),':bx','linewidth',2);
grid minor;
a = axis;
axis([1 m a(3) a(4)]);
set(gca,'FontSize',16);
legend('GPMSA',...
       'QUESO',...
       'location','north');
xlabel('Index');
ylabel('Simulation PC 2');
print -dpng simul_pc_2_both.png
waitforbuttonpress;
clf

plot(1:m,gpmsa_w(:,2)-Kmat_eta_sub0(:,2),'-bo','linewidth',2);
grid minor;
a = axis;
axis([1 m a(3) a(4)]);
set(gca,'FontSize',16);
xlabel('Index');
ylabel('Sim. PC 2, Diff');
print -dpng simul_pc_2_diff.png
waitforbuttonpress;
clf

plot(1:m,log10(abs((gpmsa_w(:,2)-Kmat_eta_sub0(:,2))./gpmsa_w(:,2))),'-bo','linewidth',2);
grid minor;
a = axis;
axis([1 m a(3) a(4)]);
set(gca,'FontSize',16);
xlabel('Index');
ylabel('Sim. PC 2, log10(|Diff/GPMSA|');
print -dpng simul_pc_2_diff_rel.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%
