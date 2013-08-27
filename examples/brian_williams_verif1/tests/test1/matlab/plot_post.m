filtChain_mh

%%%%%%%%%%%%%%%%%%%%

[f,xi] = ksdensity(sip_ip_mh_filtChain_unified(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%waitforbuttonpress;
%hold
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1');
ylabel('Marginal posterior KDE');
title('2013-08-26 - Posterior - Verif. Prob. 1');
%a=axis();
%axis([-1 1 a(3) a(4)]);

print -dpng 2013_08_26_post_param01_verif_prob_1.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%
