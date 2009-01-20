opt_data_1_;
opt_data_2_;
opt_data_3_;
%#iso_data_1_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(opt_data_1_tmpW_Time    ,opt_data_1_tmpW_Value,'--r','linewidth',2);
hold
plot(opt_data_1_initialW_Time,opt_data_1_initialW_W,'--g','linewidth',2);
plot(opt_data_1_optW_Time    ,opt_data_1_optW_W    ,'--b','linewidth',2);
%plot(iso_data_1_tmpW_Time    ,iso_data_1_tmpW_Value,'-r','linewidth',2);
%'Best optimal',...

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Evolution of mass fraction','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Reference',...
       'Initial',...
       'Optimal',...
       'location','southwest');
print -dpng opt_w.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(opt_data_1_tmpW_Time    ,opt_data_1_tmpW_Value,'--r','linewidth',2);
hold
plot(opt_data_1_initialW_Time,opt_data_1_initialW_W,'--g','linewidth',2);
plot(opt_data_1_optW_Time    ,opt_data_1_optW_W    ,'--b','linewidth',2);
a = axis;
axis([8200 10200 a(3) a(4)]);

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Evolution of mass fraction','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Reference',...
       'Initial',...
       'Optimal',...
       'location','southwest');
print -dpng opt_w_zoom.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(opt_data_1_optLambda_Time,opt_data_1_optLambda_Lambda,'-g','linewidth',2);
hold
plot(opt_data_2_optLambda_Time,opt_data_2_optLambda_Lambda,'-r','linewidth',2);
plot(opt_data_3_optLambda_Time,opt_data_3_optLambda_Lambda,'-b','linewidth',2);

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Evolution of optimal Lagrange multiplier','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Optimal Lambda 1',...
       'Optimal Lambda 2',...
       'Optimal Lambda 3',...
       'location','southwest');
print -dpng opt_lambda.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(opt_data_1_optLambda_Time,opt_data_1_optLambda_Lambda,'-g','linewidth',2);
hold
plot(opt_data_2_optLambda_Time,opt_data_2_optLambda_Lambda,'-r','linewidth',2);
plot(opt_data_3_optLambda_Time,opt_data_3_optLambda_Lambda,'-b','linewidth',2);
a = axis;
axis([8200 10200 -5 5]);

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Evolution of optimal Lagrange multiplier','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Optimal Lambda 1',...
       'Optimal Lambda 2',...
       'Optimal Lambda 3',...
       'location','southwest');
print -dpng opt_lambda_zoom.png
waitforbuttonpress
clf;
