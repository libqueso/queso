clear all;
path(path,'/h2/prudenci/DownloadedPackages/mcmc02Jun2008/code');
path(path,'/h2/prudenci/svn/pecos/uq/trunk/appls/mcmc/templateex');
uqTemplateExOutput;

figure(1); clf
mcmcplot(chainCpp,[1:4],resultsCpp.names,'chainpanel');

[nsimu,npar]=size(chainCpp);
figure(2); clf
mcmcplot(chainCpp,[1,2],resultsCpp.names,'pairs',0);
