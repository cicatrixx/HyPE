function flowdist = flowdistanceSD(acc,fdir,outbasin)
%% Loops from highest to lowest accumulation cells within the basin and 
% accumulates hydrological distance from the outlet by successively adding
% distance from downstream neighbour. Only evaluates distance for basin
% cells as loop only used inside cell index

defdirs;

% fix acc map : value outside basin is nan
acc(outbasin) = nan;

% create placeholder map
flowdist = zeros(size(acc), 'single');

% Sort acc to find order to work upstream from the outlet, i.e. decreasing acc
[~,order] = sort(acc(:), 'descend', 'MissingPlacement','first'); % nans appear first and are hence ignored by default
[nr,nc] = size(acc);

% ID basin outlet and set distance to 0 there
[~,outlet]=max(acc(:));
firstk = find(order==outlet);
%%
for k=firstk:nr*nc
    % get linear index of current point
    i = order(k); 
    
    if k==firstk % should be outlet! where dist is 0
        %assert(dir(i)==0 | dir(i)==255);
        flowdist(i) =1;
        continue;
    end
    
%      %skip if cell is outside the basin - not needed as nans ignored
%      outsid loop
%      if outbasin(i)
%          break;
%      end
%     
    % now we are on land, upstream of the outlet. Just increase the
    % distance from the downstream nbr
    r = rem(i-1, nr)+1;
    c = fix((i-1)/nr)+1;
    % find neighbour
    d = fdir(i);
    rj = min(nr,max(1,r+drow(d)));
    cj = min(nc,max(1,c+dcol(d)));
%     assert(dir(rj,cj)==0 | flowdist(rj,cj)>0);
%     if flowdist(rj,cj)==0 && dir(rj,cj)~=0
%         fdm=flowdist(rj-1:rj+1,cj-1:cj+1);
%         add=max(fdm(:));
%     else
        add=1;
%     end
    
    flowdist(r,c) = flowdist(rj,cj)+add;
end

flowdist=flowdist-1; %because the indus main outlet is set to 1 and all other outlets are also 1 and need to be set to 0. Vals set to 0 are selected by outlet generator later
%%
% cmap=jet(256);
% cmap(1,:)=[1 1 1];
% a=find(flowdist(:)==1 & acc(:)~=0);
% [r,c]=ind2sub(size(acc),a);
% figure(1);clf;imagesc(flowdist);axis image;colormap(cmap);hold on
% plot(c,r,'.r','markersize',15)
% hold off
% end