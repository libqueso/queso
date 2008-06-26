clear all;
path(path,'/h2/prudenci/DownloadedPackages/mcmc02Jun2008/code');
path(path,'/h2/prudenci/svn/pecos/uq/trunk/appls/mcmc/algaeex');
uqAlgaeExOutput;

resultsCpp.theta(1) =  0.4446;
resultsCpp.theta(2) =  0.0276;
resultsCpp.theta(3) =  0.1024;
resultsCpp.theta(4) = 11.6216;
resultsCpp.theta(5) =  0.0237;
resultsCpp.theta(6) =  1.0055;
resultsCpp.theta(7) =  3.8408;
resultsCpp.theta(8) =  1.4948;
resultsCpp.theta(9) =  9.6217;

figure
mcmcplot(chainCpp,[],resultsCpp,'pairs');
figure
mcmcplot(chainCpp,[],resultsCpp,'denspanel',2);

chainstats(chainCpp,resultsCpp);

load algaedata.mat

modelfun = @(d,th) algaefun(d(:,1),th,th(7:9),d);

% We sample 500 parameter realizations from |chain| and |s2chain|
% and calculate the predictive plots.
nsample = 500;
out = mcmcpred(resultsCpp,chainCpp,s2chainCpp,data.xdata,modelfun,nsample);
figure
mcmcpredplot(out);
% add the 'y' observations to the plot
hold on
for i=1:3
  subplot(3,1,i)
  hold on
  plot(data.ydata(:,1),data.ydata(:,i+1),'s'); 
  ylabel(''); title(data.ylabels(i+1));
  hold off
end
xlabel('days');
