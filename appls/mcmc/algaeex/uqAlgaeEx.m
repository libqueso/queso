clear all;
path(path,'/h2/prudenci/DownloadedPackages/mcmc02Jun2008/code');
path(path,'/h2/prudenci/svn/pecos/uq/trunk/appls/mcmc/algaeex');
uqAlgaeExOutput;

figure
mcmcplot(chainCpp,[],resultsCpp,'pairs');
figure
mcmcplot(chainCpp,[],resultsCpp,'denspanel',2);

chainstats(chainCpp,resultsCpp);

