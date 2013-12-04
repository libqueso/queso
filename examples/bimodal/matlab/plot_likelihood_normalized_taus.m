%cd ../test/test_2013_12_02/outputData
cd outputData

tau1 = 0.00585938; %in display_sub0.txt, look for first occurence of prevExponent after 'generateSequence(), level 2, step 3'
tau2 = 0.0524597;  %in display_sub0.txt, look for first occurence of prevExponent after 'generateSequence(), level 3, step 3'
tau3 = 0.348566;   %in display_sub0.txt, look for first occurence of prevExponent after 'generateSequence(), level 4, step 3'
tau4 = 1;

x=-50.:.1:200.;
[m,n]=size(x);
deltaX = (max(x)-min(x))/(n-1);
y1=(x-10).*(x-10)/(2*1);
y2=(x-100).*(x-100)/(2*25);
z1=(exp(-y1))/(1*sqrt(2*pi));
z2=(exp(-y2))/(5*sqrt(2*pi));

% Level 1------------------------------------------------------------------
z=((z1+ z2)/2).^tau1;
intZ = 0.;
for i=1:(n-1)
intZ = intZ + z(i)*deltaX;
end
intZ=intZ;

fprintf(1,'Level 1, area before normalization = %4.3f,',intZ);
z = z/intZ;
intZ = 0.;
for i=1:(n-1)
intZ = intZ + z(i)*deltaX;
end
intZ=intZ;
fprintf(1,'\t area after normalization = %4.3f\n',intZ);

x1=x;
f1=z;
clear z; clear intZ; 

% Level 2------------------------------------------------------------------
z=((z1+ z2)/2).^tau2;
intZ = 0.;
for i=1:(n-1)
intZ = intZ + z(i)*deltaX;
end
intZ=intZ;

fprintf(1,'Level 2, area before normalization = %4.3f,',intZ);
z = z/intZ;
intZ = 0.;
for i=1:(n-1)
intZ = intZ + z(i)*deltaX;
end
intZ=intZ;
fprintf(1,'\t area after normalization = %4.3f\n',intZ);

x2=x;
f2=z;
clear z; clear intZ; 

% Level 3------------------------------------------------------------------
z=((z1+ z2)/2).^tau3;
intZ = 0.;
for i=1:(n-1)
intZ = intZ + z(i)*deltaX;
end
intZ=intZ;

fprintf(1,'Level 3, area before normalization = %4.3f,',intZ);
z = z/intZ;
intZ = 0.;
for i=1:(n-1)
intZ = intZ + z(i)*deltaX;
end
intZ=intZ;
fprintf(1,'\t area after normalization = %4.3f\n',intZ);

x3=x;
f3=z;
clear z; clear intZ;

% Level 4------------------------------------------------------------------
z=((z1+ z2)/2).^tau4;
intZ = 0.;
for i=1:(n-1)
intZ = intZ + z(i)*deltaX;
end
intZ=intZ;

fprintf(1,'Level 4, area before normalization = %4.3f,',intZ);
z = z/intZ;
intZ = 0.;
for i=1:(n-1)
intZ = intZ + z(i)*deltaX;
end
intZ=intZ;
fprintf(1,'\t area after normalization = %4.3f\n',intZ);

x4=x;
f4=z;
% -------------------------------------------------------------------------
plot(x1,f1, '--r',x2,f2, '--g',x3,f3,'--b',x4,f4,'--c','linewidth',3)
set(gca,'fontsize',20);
a=axis;
axis([-50 200 a(3) .15]);
legend(['\tau = ', num2str(tau1)],...
       ['\tau = ', num2str(tau2)],...
       ['\tau = ', num2str(tau3)],...
       ['\tau = ', num2str(tau4)],...
        'location','northeast');
ylabel('f^{\tau}(\theta)','fontsize',16);
xlabel('\theta','fontsize',16);
title('Intermediary normalized likelihood functions','fontname', 'Times', 'fontsize',20);
grid minor
set(gca,'fontsize',16);

print -dpng bimodal_likelihood_taus_normalized.png

%cd ../../../matlab
cd ..
