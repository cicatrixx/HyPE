function [ro,co,sites] = find_outletsInChannel(channels,flowdist,do)
% do = number of cells in between outlets
% Find outlets based on hydrological distance map
sites = find(channels & rem(flowdist, do)==0 ); % select pixels at equal distance intervals AND within channel network AND WITHIN BASIN BOUNDS
[ro,co] = ind2sub(size(flowdist), sites);

% figure(1);imagesc(X); hold on; 
% plot(co,ro,'w.', 'markersize',20); 
% hold off;
end