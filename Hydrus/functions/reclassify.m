
%% Reclassify data based on breakpoints provided
function idata_reclass= reclassify(idata, breakpoints,newclass, classname)
% Start and end points are auto generated so breakpoints only need to
% include inner bounding value. min and max bounds evaluated based on the data
% Reclassification done assuming lbounds(i) <= idata < % ubounds(i) 
% except for last class where ubounds is not applied
% Plot is made if classnames is defined
% Nans are kept nans. 0s might be reclassed as class1 so be careful
idata_reclass= 0*idata; 

breakpoints= reshape(breakpoints,[numel(breakpoints),1]); % read breakpoints as column vector
lbounds=[min(idata(:)); breakpoints];
ubounds= [breakpoints; max(idata(:))];

for i=1:(length(lbounds)-1 )
    idata_reclass(idata>=lbounds(i)& idata<ubounds(i))=newclass(i); 
end

% For last class, no ubound
i=i+1;
idata_reclass(idata>=lbounds(i))=newclass(i); 

if exist('classname','var')
    figure;
    subplot(1,2,1)
    imagescnan(idata_reclass)
    subplot(1,2,2)
    histogram(idata_reclass(idata_reclass>0))
    grid on
    xticklabels(classname)
end