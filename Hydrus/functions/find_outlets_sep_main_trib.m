function [ro,co,channels] = find_outlets_sep_main_trib(flowdist,acc,do)

%allIndus=load(fullfile(pwd,'data/ASIA/Basin/Basin_5.mat'), 'acc');
%ma = mean(allIndus.acc(:)); %to standardize channel definition

% Find outlets based on hydrological distance from main outlet
ma = mean(acc(:)); %to standardize channel definition
channels = acc>ma; % only cells above mean flow accumulation area

stations = find(channels & rem(flowdist, do)==0 & flowdist>0); % select pixels at equal distance intervals AND within channel network AND WITHIN BASIN BOUNDS
[ro,co] = ind2sub(size(acc), stations);

% figure(1);imagesc(X); hold on; 
% plot(co,ro,'w.', 'markersize',20); 
% hold off;
end
