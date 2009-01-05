%string1 = sprintf('Refer: A=%9.3e, E=%9.3e',case5_refA,case6_refE);
%string2 = sprintf('Guess: A=%9.3e, E=%9.3e',case5_guessA,case6_guessE);

case5_;
case6_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case5_refW_Time,case5_refW_Value,'+k','linewidth',2);
hold
plot(case5_w_Time,case5_w_Value,'-b','linewidth',2);
plot(case6_w_Time,case6_w_Value,'-r','linewidth',2);

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Evolution of mass fraction','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Reference',...
       'Case5',...
       'Case6',...
       'location','southwest');
print -dpng case5_case6_w.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case5_refW_Time,case5_refW_Value,'+k','linewidth',3);
hold
plot(case5_w_Time,case5_w_Value, '-b','linewidth',3);
plot(case6_w_Time,case6_w_Value, '-r','linewidth',3);
a = axis;
axis([8000 9600 a(3) a(4)]);

ylabel('w','fontsize',20);
xlabel('time (s)','fontsize',20);
title('Zoom of the evolution of mass fraction','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Reference',...
       'Case5',...
       'Case6',...
       'location','southwest');
print -dpng case5_case6_w_zoom.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case5_wA_Time,case5_wA_Value,'-b','linewidth',3);
hold
plot(case6_wA_Time,case6_wA_Value,'-r','linewidth',3);

ylabel('w_A','fontsize',20);
xlabel('time (s)','fontsize',20);
title('w_A','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Case5',...
       'Case6',...
       'location','southwest');
print -dpng case5_case6_wA.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case5_wE_Time,case5_wE_Value,'-b','linewidth',3);
hold
plot(case6_wE_Time,case6_wE_Value,'-r','linewidth',3);

ylabel('w_E','fontsize',20);
xlabel('time (s)','fontsize',20);
title('w_E','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Case5',...
       'Case6',...
       'location','southwest');
print -dpng case5_case6_wE.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case5_lambda_Time,case5_lambda_Value,'-b','linewidth',3);
hold
plot(case6_lambda_Time,case6_lambda_Value,'-r','linewidth',3);

ylabel('\lambda','fontsize',20);
xlabel('time (s)','fontsize',20);
title('\lambda','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Case5',...
       'Case6',...
       'location','southwest');
print -dpng case5_case6_lambda.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case5_lambdaA_Time,case5_lambdaA_Value,'-b','linewidth',3);
hold
plot(case6_lambdaA_Time,case6_lambdaA_Value,'-r','linewidth',3);

ylabel('\lambda_A','fontsize',20);
xlabel('time (s)','fontsize',20);
title('\lambda_A','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Case5',...
       'Case6',...
       'location','southwest');
print -dpng case5_case6_lambdaA.png
waitforbuttonpress
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(case5_lambdaE_Time,case5_lambdaE_Value,'-b','linewidth',3);
hold
plot(case6_lambdaE_Time,case6_lambdaE_Value,'-r','linewidth',3);

ylabel('\lambda_E','fontsize',20);
xlabel('time (s)','fontsize',20);
title('\lambda_E','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('Case5',...
       'Case6',...
       'location','southwest');
print -dpng case5_case6_lambdaE.png
waitforbuttonpress
clf;

