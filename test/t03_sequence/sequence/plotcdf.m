load sortedData.dat;
load cdf.dat;
load cdfup.dat;
load cdfdown.dat;
[aa bb]=size(sortedData)
for jj=1:bb
    figure(jj);plot(sortedData(:,jj),cdf(:,jj),sortedData(:,jj),cdfup(:,jj),sortedData(:,jj),cdfdown(:,jj));
end
