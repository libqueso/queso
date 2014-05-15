%cd ../test/test_2013_12_02/outputData
cd outputData

tau1 = 0.00585938; %in display_sub0.txt, look for first occurence of prevExponent after 'generateSequence(), level 2, step 3'
tau2 = 0.0524597;  %in display_sub0.txt, look for first occurence of prevExponent after 'generateSequence(), level 3, step 3'
tau3 = 0.348566;   %in display_sub0.txt, look for first occurence of prevExponent after 'generateSequence(), level 4, step 3'
tau4 = 1;


x=-100.:.1:300.;
y1=(x-10).*(x-10)/(2*1);
y2=(x-100).*(x-100)/(2*25);
z1=(exp(-y1))/(1*sqrt(2*pi));
z2=(exp(-y2))/(5*sqrt(2*pi));
y=(z1+z2)/2.;

plot(x,y.^tau1,'r-','linewidth',2);
hold
plot(x,y.^tau2,'g-','linewidth',2);
plot(x,y.^tau3,'b-','linewidth',2);
plot(x,y.^tau4,'c-','linewidth',2);

a=axis;
axis([-100 300 a(3) 1]);
ylabel('f^{\tau}(\theta)','fontsize',16);
xlabel('\theta','fontsize',16);

title('Intermediary likelihood functions','fontname', 'Times', 'fontsize',20);
grid minor
set(gca,'fontsize',16);

legend(['\tau = ', num2str(tau1)],...
       ['\tau = ', num2str(tau2)],...
       ['\tau = ', num2str(tau3)],...
       ['\tau = ', num2str(tau4)],...
        'location','northeast');

print -dpng bimodal_likelihood_taus.png
hold off;

%cd ../../../matlab
cd ..
