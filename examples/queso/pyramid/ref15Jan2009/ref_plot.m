ref_data_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(ref_data_tmpW1_Time,ref_data_tmpW1_W   ,'-b' ,'linewidth',2);
hold
plot(ref_data_tmpW2_Time,ref_data_tmpW2_W   ,'--b','linewidth',2);
plot(ref_data_tmpW3_Time,ref_data_tmpW3_W   ,'-r' ,'linewidth',2);
plot(ref_data_tmpW4_Time,ref_data_tmpW4_W   ,'--r','linewidth',2);
plot(ref_data_tmpW_Time ,ref_data_tmpW_Value,'--g','linewidth',2);

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Evolution of mass fraction','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Material 1',...
       'Material 2',...
       'Material 3',...
       'Material 4',...
       'Composed Sample',...
       'location','southwest');
print -dpng ref_w.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(ref_data_tmpW1_Time,ref_data_tmpW1_W   ,'-b' ,'linewidth',2);
hold
plot(ref_data_tmpW2_Time,ref_data_tmpW2_W   ,'--b','linewidth',2);
plot(ref_data_tmpW3_Time,ref_data_tmpW3_W   ,'-r' ,'linewidth',2);
plot(ref_data_tmpW4_Time,ref_data_tmpW4_W   ,'--r','linewidth',2);
plot(ref_data_tmpW_Time ,ref_data_tmpW_Value,'--g','linewidth',2);
a = axis;
axis([8200 10200 a(3) a(4)]);

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Evolution of mass fraction','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Material 1',...
       'Material 2',...
       'Material 3',...
       'Material 4',...
       'Composed Sample',...
       'location','southwest');
print -dpng ref_w_zoom.png
waitforbuttonpress
clf;

