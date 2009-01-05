%string1 = sprintf('Refer: A=%9.3e, E=%9.3e',case1_refA,case1_refE);
%string2 = sprintf('Guess: A=%9.3e, E=%9.3e',case1_guessA,case1_guessE);

case1_;
case2_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case1_refW_Time,case1_refW_Value,'--k','linewidth',2);
hold
plot(case1_w_Time,case1_w_Value,'-b','linewidth',2);
plot(case2_w_Time,case2_w_Value,'-r','linewidth',2);

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Evolution of mass fraction','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Reference',...
       'Case1',...
       'Case2',...
       'location','southwest');
print -dpng case1_case2_w.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case1_refW_Time,case1_refW_Value,'--k','linewidth',3);
hold
plot(case1_w_Time,case1_w_Value, '-b','linewidth',3);
plot(case2_w_Time,case2_w_Value, '-r','linewidth',3);
a = axis;
axis([8000 9600 a(3) a(4)]);

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Zoom of the evolution of mass fraction','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Reference',...
       'Case1',...
       'Case2',...
       'location','southwest');
print -dpng case1_case2_w_zoom.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case1_wA_Time,case1_wA_Value,'-b','linewidth',3);
hold
plot(case2_wA_Time,case2_wA_Value,'-r','linewidth',3);

ylabel('w_A','fontsize',20);
xlabel('time (s)','fontsize',20);
title('w_A','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Case1',...
       'Case2',...
       'location','southwest');
print -dpng case1_case2_wA.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case1_wE_Time,case1_wE_Value,'-b','linewidth',3);
hold
plot(case2_wE_Time,case2_wE_Value,'-r','linewidth',3);

ylabel('w_E','fontsize',20);
xlabel('time (s)','fontsize',20);
title('w_E','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Case1',...
       'Case2',...
       'location','southwest');
print -dpng case1_case2_wE.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case1_lambda_Time,case1_lambda_Value,'-b','linewidth',3);
hold
plot(case2_lambda_Time,case2_lambda_Value,'-r','linewidth',3);

ylabel('\lambda','fontsize',20);
xlabel('time (s)','fontsize',20);
title('\lambda','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Case1',...
       'Case2',...
       'location','southwest');
print -dpng case1_case2_lambda.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case1_lambdaA_Time,case1_lambdaA_Value,'-b','linewidth',3);
hold
plot(case2_lambdaA_Time,case2_lambdaA_Value,'-r','linewidth',3);

ylabel('\lambda_A','fontsize',20);
xlabel('time (s)','fontsize',20);
title('\lambda_A','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Case1',...
       'Case2',...
       'location','southwest');
print -dpng case1_case2_lambdaA.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case1_lambdaE_Time,case1_lambdaE_Value,'-b','linewidth',3);
hold
plot(case2_lambdaE_Time,case2_lambdaE_Value,'-r','linewidth',3);

ylabel('\lambda_E','fontsize',20);
xlabel('time (s)','fontsize',20);
title('\lambda_E','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Case1',...
       'Case2',...
       'location','southwest');
print -dpng case1_case2_lambdaE.png
waitforbuttonpress
clf;

