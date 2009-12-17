load xdata.dat;
load cdf.dat;
load cdfup.dat;
load cdfdown.dat;
[aa bb]=size(xdata);
for jj=1:bb
    figure(jj);plot(xdata(:,jj),cdf(:,jj),xdata(:,jj),cdfup(:,jj),xdata(:,jj),cdfdown(:,jj));
end