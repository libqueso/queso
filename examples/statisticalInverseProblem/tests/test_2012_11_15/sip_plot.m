%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd outputData
sipOutput_sub0
file_sip_raw
appl_output_sub0

[nsimu,temp]=size(sip_appl_d_sub0);
[temp,npar]=size(ip_mh_rawChain_unified);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot samples of first two parameters on the plane
%
% Acknowledgments:
% ===============
% Figure 'parameters_samples_plane.png' has the purpose 
% of being compared with the output of the normal test 
% at the MCMC toolbox for MATLAB, available at 
% www.helsinki.fi/~mjlaine/mcmc/.
% Some of the commands below indeed resamble part of
% the code available at such toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c50  = chi2inv(0.50,npar);
c95  = chi2inv(0.95,npar);
cc50 = sum(sip_appl_d_sub0<c50)./nsimu;
cc95 = sum(sip_appl_d_sub0<c95)./nsimu;

plot(ip_mh_rawChain_unified(:,1),ip_mh_rawChain_unified(:,2),'.');
set(gca,'fontsize',20);
xlabel('\theta_1');
ylabel('\theta_2');
title(sprintf('Rejected = %.1f%%, c50 = %.1f%%, c95 = %.1f%%', ...
              ip_mh_rejected*100, cc50*100, cc95*100))
hold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 50% and 95% ellipses on the same plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c50  = chi2inv(0.50,2);
c95  = chi2inv(0.95,2);

t = linspace(0,2*pi)'; %'
R = chol(c50*sip_appl_covMatrix_sub0(1:2,1:2));
x = sip_appl_paramMeans_sub0(1) + R(1,1).*cos(t);
y = sip_appl_paramMeans_sub0(2) + R(1,2).*cos(t) + R(2,2).*sin(t);
plot(x,y,'r--','LineWidth',2);

R = chol(c95*sip_appl_covMatrix_sub0(1:2,1:2));
x = sip_appl_paramMeans_sub0(1) + R(1,1).*cos(t);
y = sip_appl_paramMeans_sub0(2) + R(1,2).*cos(t) + R(2,2).*sin(t);
plot(x,y,'r--','LineWidth',2);
grid minor;
print -dpng parameters_samples_plane.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot marginal posterior pdfs of all 4 parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KDE plots ---------------------------------------------------------------
[f1,x1] = ksdensity(ip_mh_rawChain_unified(1,:),'function','pdf');
plot(x1,f1,'-b','linewidth',3)
hold;
[f2,x2] = ksdensity(ip_mh_rawChain_unified(2,:),'function','pdf');
plot(x2,f2,'-r','linewidth',3)
[f3,x3] = ksdensity(ip_mh_rawChain_unified(3,:),'function','pdf');
plot(x3,f3,'--b','linewidth',3)
[f4,x4] = ksdensity(ip_mh_rawChain_unified(4,:),'function','pdf');
plot(x4,f4,'--r','linewidth',3)
hold;
ylabel('marginal pdf','fontsize',20);
xlabel('Parameter values','fontsize',20);
title('Posterior marginal pdfs of parameters','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('param1',...
       'param2',...
       'param3',...
       'param4',...
       'location','northwest');
print -dpng parameters_PDF.png

cd ..
