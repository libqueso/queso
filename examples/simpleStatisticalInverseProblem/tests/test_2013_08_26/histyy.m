function varargout = histyy(data1,nBins1,data2,nBins2,varargin)
%% function to plot two histograms on the same plot on two different axes.
% Syntax
% ------
% histyy(data1, nBins1, data2, nBins2)
% histyy(data1, nBins1, data2, nBins2)
% histyy(data1, nBins1, data2, nBins2, label1, label2)
% hAx = histyy(...)
%
% Inputs
% ------
% data1 - vector containing first histogram
% nBins1 - number of bins in first histogram
% data2 - vector containing second histogram
% nBins2 - number of bins in second histogram
% label1 - ylabel for first histogram
% label2 - ylabel for second histogram
%
% Outputs
% -------
% hAx - vector containing handles to axes
% hAx(1) contains handle to first histogram axis
% hAx(2) contains handle to second histogram axis
%
% Author : araja 11/07/11

if isempty(varargin)
label1 = '';
label2 = '';
elseif length(varargin) ~= 2
error('Provide both label strings.');
else
label1 = varargin{1};
label2 = varargin{2};
end

if isempty(nBins1)
nBins1 = 10;
end
if isempty(nBins2)
nBins2 = 10;
end

% create one axis or get current axes
hAx1 = gca;

% get position of axis
posAx1 = get(hAx1, 'Position');

% create an overlapping axis at the same location
hAx2 = axes('Position', posAx1);

% histogram for first data vector
hist(hAx1,data1,nBins1);
set(findobj(hAx1,'Type','patch'),'FaceColor','b','EdgeColor','w', 'facealpha',0.75);


% histogram for second data vector
hist(hAx2, data2,nBins2);

% make second axis transparent
set(hAx2,'Color','none');

% change color of second histogram
set(findobj(hAx2,'Type','patch'),'FaceColor','g','EdgeColor','w', 'facealpha',0.75);

% ylabel for histogram 1
ylabel(hAx1,label1);

% ylabel for histogram 2
set(hAx2,'YAxisLocation','right');
ylabel(hAx2,label2);

% set x-axis limits
lim1 = get(hAx1, 'XLim');
lim2 = get(hAx2, 'XLim');

set(hAx1, 'XLim', [min([lim1(1) lim2(1)]) max([lim1(2) lim2(2)])]);
set(hAx2, 'XLim', [min([lim1(1) lim2(1)]) max([lim1(2) lim2(2)])]);
set(hAx2, 'XTickLabel',[])

% output
if nargout == 1
varargout{1} = [hAx1 hAx2];
end


