cd outputData
sipOutput_sub0
sfpOutput_sub0

plot(ip_mh_rawChain_corrViaFftLags_sub0,ip_mh_rawChain_corrViaFftInitPos0_sub0(1,:),'-b','linewidth',2);
hold
plot(ip_mh_filtChain_corrViaFftLags_sub0,ip_mh_filtChain_corrViaFftInitPos0_sub0(1,:),'-r','linewidth',2);
ylabel('Autocorrelation for \theta_1','fontsize',20);
xlabel('Lag','fontsize',20);
%title('Fig 1, Autocorrelation for \theta_1 samples','fontsize',20);
a = axis;
axis([a(1) a(2) -0.1 1]);
grid minor;
set(gca,'fontsize',20);
legend('raw chain',...
       'filtered chain',...
       'location','northeast');
print -dpng paper_plot1.png
waitforbuttonpress;
clf;

plot(ip_mh_rawChain_corrViaFftLags_sub0,ip_mh_rawChain_corrViaFftInitPos0_sub0(2,:),'-b','linewidth',2);
hold
plot(ip_mh_filtChain_corrViaFftLags_sub0,ip_mh_filtChain_corrViaFftInitPos0_sub0(2,:),'-r','linewidth',2);
ylabel('Autocorrelation for \theta_2','fontsize',20);
xlabel('Lag','fontsize',20);
%title('Fig 2, Autocorrelation for \theta_2 samples','fontsize',20);
a = axis;
axis([a(1) a(2) -0.1 1]);
grid minor;
set(gca,'fontsize',20);
legend('raw chain',...
       'filtered chain',...
       'location','northeast');
print -dpng paper_plot2.png
waitforbuttonpress;
clf;

plot(ip_mh_filtChain_unifGkdePosits_sub0(1,:),ip_mh_filtChain_unifGkdeValues_sub0(1,:),'-b','linewidth',2);
hold
x = ip_mh_filtChain_unifGkdePosits_sub0(1,:);
plot(x,(exp(-(x+1).*(x+1)/8))/2/sqrt(2*pi),'--r','linewidth',2);
ylabel('Posterior marginal pdf','fontsize',20);
xlabel('\theta_1','fontsize',20);
%title('Fig 3, Pdfs for \theta_1','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('QUESO',...
       'Analytic',...
       'location','northwest');
print -dpng paper_plot3.png
waitforbuttonpress;
clf;

plot(ip_mh_filtChain_unifGkdePosits_sub0(2,:),ip_mh_filtChain_unifGkdeValues_sub0(2,:),'-b','linewidth',2);
hold
x = ip_mh_filtChain_unifGkdePosits_sub0(2,:);
plot(x,(exp(-(x-2).*(x-2)/2))/sqrt(2*pi),'--r','linewidth',2);
ylabel('Posterior marginal pdf','fontsize',20);
xlabel('\theta_2','fontsize',20);
%title('Fig 4, Pdfs for \theta_2','fontsize',20);
grid minor;
set(gca,'fontsize',20);
legend('QUESO',...
       'Analytic',...
       'location','northwest');
print -dpng paper_plot4.png
waitforbuttonpress;
clf;

plot(fp_mc_QoiSeq_unifGkdePosits_sub0(1,:),fp_mc_QoiSeq_unifGkdeValues_sub0(1,:),'-b','linewidth',2);
ylabel('Pdf','fontsize',20);
xlabel('QoI = \theta_1 + \theta_2','fontsize',20);
%title('Fig 5, Pdf for QoI','fontsize',20);
a = axis;
axis([-9 11 a(3) a(4)]);
grid minor;
set(gca,'fontsize',20);
%legend('QUESO',...
%       'Analytic',...
%       'location','northwest');
print -dpng paper_plot5.png
