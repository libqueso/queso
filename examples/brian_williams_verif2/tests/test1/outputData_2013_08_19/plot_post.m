filtChain_mh

muVec = [2.24522
5.61304
]; 
sigmaMat = [862.073 -344.817
-344.817 137.957
];

%%%%%%%%%%%%%%%%%%%%

[f,xi] = ksdensity(sip_ip_mh_filtChain_unified(:,1),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%waitforbuttonpress;
%hold
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1');
ylabel('Marginal posterior KDE');
title('2013-08-19 - Posterior - Verif. Prob. 2');
%a=axis();
%axis([-1 1 a(3) a(4)]);

print -dpng 2013_08_19_post_param01_verif_prob_2.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

[f,xi] = ksdensity(sip_ip_mh_filtChain_unified(:,2),'function','pdf');
plot(xi,f,'-b','linewidth',2)
%waitforbuttonpress;
%hold
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 2');
ylabel('Marginal posterior KDE');
title('2013-08-19 - Posterior - Verif. Prob. 2');
%a=axis();
%axis([0 1000 a(3) a(4)]);

print -dpng 2013_08_19_post_param02_verif_prob_2.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%

path(path,'../matlab')

data=[sip_ip_mh_filtChain_unified(:,1), sip_ip_mh_filtChain_unified(:,2)];
[bandwidth,density,X,Y]=kde2d(data,64);
surf(X,Y,density);
%waitforbuttonpress;
%hold
grid minor;
set(gca,'FontSize',16);
xlabel('Parameter 1 ()');
ylabel('Parameter 2 ()');
zlabel('Marginal 2D posterior KDE');
%a=axis();
%axis([5 40 4 22 a(5) a(6)]);

view([315,45]);
title('2013-08-19 - Posterior - Verif. Prob. 2 - View A');
print -dpng 2013_08_19_post_kde_params_1_and_2_verif_prob_2_viewA.png
waitforbuttonpress;
view([225,45]);
title('2013-08-19 - Posterior - Verif. Prob. 2 - View B');
print -dpng 2013_08_19_post_kde_params_1_and_2_verif_prob_2_viewB.png
waitforbuttonpress;
view([135,45]);
title('2013-08-19 - Posterior - Verif. Prob. 2 - View C');
print -dpng 2013_08_19_post_kde_params_1_and_2_verif_prob_2_viewC.png
waitforbuttonpress;
view([45,45]);
title('2013-08-19 - Posterior - Verif. Prob. 2 - View D');
print -dpng 2013_08_19_post_kde_params_1_and_2_verif_prob_2_viewD.png
waitforbuttonpress;
view([0,90]);
title('2013-08-19 - Posterior - Verif. Prob. 2 - View E');
print -dpng 2013_08_19_post_kde_params_1_and_2_verif_prob_2_viewE.png
waitforbuttonpress;
clf;

%%%%%%%%%%%%%%%%%%%%
