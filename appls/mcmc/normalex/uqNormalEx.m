clear all;
path(path,'/h2/prudenci/DownloadedPackages/mcmc02Jun2008/code');
path(path,'/h2/prudenci/svn/pecos/uq/trunk/appls/mcmc/normalex');
uqNormalExOutput;

figure(1); clf
mcmcplot(chainCpp,[1:4],resultsCpp.names,'chainpanel');

[nsimu,npar]=size(chainCpp);
c50 = chiqf(0.50,npar);
c95 = chiqf(0.95,npar);
cc50 = sum(dCpp<c50)./nsimu;
cc95 = sum(dCpp<c95)./nsimu;
figure(2); clf
mcmcplot(chainCpp,[1,2],resultsCpp.names,'pairs',0);

title(sprintf('Rejected = %.1f%%, c50 = %.1f%%, c95 = %.1f%%', ...
              resultsCpp.rejected*100, cc50*100, cc95*100));

c50  = chiqf(0.50,2);
c95  = chiqf(0.95,2);

hold on
% The code below assumes that the mahalanobisMatrixCpp,
% used for the computations of malananobis distances,
% is the covariance matrix.
ellipse(priorMeanValuesCpp(1:2),c50*mahalanobisMatrixCpp(1:2,1:2),'r--','LineWidth',2);
ellipse(priorMeanValuesCpp(1:2),c95*mahalanobisMatrixCpp(1:2,1:2),'r-','LineWidth',2);
axis equal
