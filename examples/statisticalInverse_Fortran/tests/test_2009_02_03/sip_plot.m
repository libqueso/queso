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
% Figure 'fig1.png' has the purpose of being compared with
% the output of the normal test at the MCMC toolbox for
% MATLAB, available at www.helsinki.fi/~mjlaine/mcmc/.
% Some of the commands below indeed resamble part of
% the code available at such toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c50  = chi2inv(0.50,npar);
c95  = chi2inv(0.95,npar);
cc50 = sum(sip_appl_d_sub0<c50)./nsimu;
cc95 = sum(sip_appl_d_sub0<c95)./nsimu;

plot(ip_mh_rawChain_unified(:,1),ip_mh_rawChain_unified(:,2),'.');
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

print -dpng fig1.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot marginal posterior pdfs of all 4 parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(ip_mh_rawChain_unifGkdePosits_sub0(1,:),ip_mh_rawChain_unifGkdeValues_sub0(1,:),'-b');
hold
plot(ip_mh_rawChain_unifGkdePosits_sub0(2,:),ip_mh_rawChain_unifGkdeValues_sub0(2,:),'-r');
plot(ip_mh_rawChain_unifGkdePosits_sub0(3,:),ip_mh_rawChain_unifGkdeValues_sub0(3,:),'--b');
plot(ip_mh_rawChain_unifGkdePosits_sub0(4,:),ip_mh_rawChain_unifGkdeValues_sub0(4,:),'--r');

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
print -dpng fig2.png
