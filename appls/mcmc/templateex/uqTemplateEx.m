clear all;
path(path,getenv('MCMC_TOOLBOX_PATH_FOR_PECOS_TOOLKIT'));
path(path,pwd);
uqTemplateExOutput;

figure(1); clf
mcmcplot(chainCpp,[1:4],resultsCpp.names,'chainpanel');

[nsimu,npar]=size(chainCpp);
figure(2); clf
mcmcplot(chainCpp,[1,2],resultsCpp.names,'pairs',0);
