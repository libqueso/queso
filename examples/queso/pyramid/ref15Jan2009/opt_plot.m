opt_data_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(opt_data_tmpW_Time,opt_data_tmpW_Value,'--g','linewidth',2);
hold
plot(opt_data_optW_Time,opt_data_optW_W    ,'--b','linewidth',2);

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Evolution of mass fraction','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Referece',...
       'Optimal',...
       'location','southwest');
print -dpng opt_w.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(opt_data_tmpW_Time,opt_data_tmpW_Value,'--g','linewidth',2);
hold
plot(opt_data_optW_Time,opt_data_optW_W    ,'--b','linewidth',2);
a = axis;
axis([8200 10200 a(3) a(4)]);

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Evolution of mass fraction','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Referece',...
       'Optimal',...
       'location','southwest');
print -dpng opt_w_zoom.png
waitforbuttonpress
clf;

