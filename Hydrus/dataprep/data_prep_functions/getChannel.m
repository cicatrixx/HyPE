function [rch,cch,ich,channel]=getChannel(acc, outside)
% Returns row, column and index number of cells that are channel>nanmean(acc(:))

% set not basin cells to nan
acc(outside)=nan;

% Find channel
ma = nanmean(acc(:));
channel = acc>ma & ~outside; % only cells above mean flow accumulation area and inside basin
ich=find(channel);
[rch,cch]=ind2sub(size(acc), ich);

end