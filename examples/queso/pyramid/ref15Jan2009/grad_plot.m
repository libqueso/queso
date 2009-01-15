grad_data_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(grad_data_tmpW_Time,grad_data_tmpW_Value,'--g','linewidth',2);
hold
plot(grad_data_optW_Time,grad_data_optW_W    ,'--b','linewidth',2);

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Evolution of mass fraction','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Referece',...
       'Optimal',...
       'location','southwest');
print -dpng grad_w.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(grad_data_tmpW_Time,grad_data_tmpW_Value,'--g','linewidth',2);
hold
plot(grad_data_optW_Time,grad_data_optW_W    ,'--b','linewidth',2);
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
print -dpng grad_w_zoom.png
waitforbuttonpress
clf;

