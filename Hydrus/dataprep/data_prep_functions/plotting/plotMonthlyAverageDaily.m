function [avg]= plotMonthlyAverageDaily(data,skip,createplot,figtitles,ylabs)
% Reshapes matrix w/ TS values along columns to plot all yearly data as
% overlapping lines and a bold line indicating monthly averages for all
% period. Each row is plotted as a separate subplot
%   data: monthly values in columns and different datasets in rows
%   skip: # years from the start to skip
%   plotnew: create plot or not
%   figtitles, ylab: titles and units of Y data. can be single val or a vector

%   avg:  monthly avg values: ndata x nmonths

nRow = size(data,1);
avg = zeros(nRow,12);
cmap = hsv(nRow);
if createplot ==1;figure;end
for r= 1:nRow
    data12=[];
    data12 = reshape(data(r,skip*12+1:end),12,size(data,2)/12); %months X years
    avg(r,:) = nanmean(data12(:,skip+1:end),2)';
    if createplot ==1
        subplot(ceil(nRow/2),2,r)
        hold all
        l1=plot(data12,'LineWidth',.5,'Color',[0.5 0.5 0.5]);
        l2=plot(avg(r,:),'r','LineWidth',2,'Color',cmap(r,:));
        if length(figtitles)>1
            ylabel(ylabs(r))        
            title(figtitles{r})
        elseif length(figtitles)==1
            ylabel(ylabs)
            title(figtitles)
        end
        %set(gca,'XTick',1:1:12,'XTickLabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
        set(gca,'XTick',1:1:12,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'})
        xlim([1 12])
        box on
        grid on
        legend([l1(1),l2],{'For each year' ,'Long-term Avg'})
    end
end
end