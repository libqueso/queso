%rawChain_ml
nada

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,1),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 1');
print -dpng param_01.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,2),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 2');
print -dpng param_02.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,3),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 3');
print -dpng param_03.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,4),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 4');
print -dpng param_04.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,5),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 5');
print -dpng param_05.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,6),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 6');
print -dpng param_06.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,7),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 7');
print -dpng param_07.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,8),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 8');
print -dpng param_08.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,9),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 9');
print -dpng param_09.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,10),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 10');
print -dpng param_10.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,11),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 11');
print -dpng param_11.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,12),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 12');
print -dpng param_12.png
waitforbuttonpress;
clf

%%%%%%%%%%%%%%%%%%%%%%

[g,gi] = ksdensity(gcm_ml_26_rawChain_unified(32000:100:132000,13),'function','pdf');
plot(gi,g,'-r','linewidth',2);
title('Parameter 13');
print -dpng param_13.png
waitforbuttonpress;
clf


