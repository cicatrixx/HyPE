function basin = fastfindupstream_DisSD4(acc,fdir,adir,flowdist,drow,dcol,outlet,showplot, dowin)
%% fast-find catchment upstream of given outlet position upto dowin cells upstream of cell
% Updated on 23 Mar 2020 to use findmainstream to get dowin automatically
% Can use masked acc because of sort descend
% Changes made
    % corrected indexing
    % skip cell downstream as well as outside search radius
    % faster init of basin
    % calculate distance threshold as furthest point in basin
    % take min on calculated dist or 2500 cells as dowin
    % add evalcell to check what cells are evaluated
    % add showplot
    % filter order based on fdir and flowdist prior to loop
    % SELECTED FOR TESTING W HYDRUS03
 
StartDis = flowdist(outlet);% Start position
if ~exist('dowin','var')    
    strm = find_mainstream(acc,outlet,adir);
    [~,seli] = min(acc(strm)); % idx in mainstrem vector
    mywin = flowdist(strm(seli))- StartDis;
    dowin= min(2500, mywin);
end

[nr,nc] = size(acc);
[~,order_all] = sort(acc(:),'descend');

% Select only cells within search radius
% Dis_order = Dis(order);

basin = zeros(nr,nc, 'logical'); % reset at empty map
evalcell = zeros(nr,nc, 'logical'); % tracker to check if cell is evaluated
basin(outlet) = true; % outlet is lowest grid cell in catchment

% Select only cells that lie within the full basin, but upstream and within
% dowin of the outlet
% fdir(1748762)=4; % temporarily set basin downstream outlet fdir to 4 as
% it is also =9! Not necessary cos it is ok if basin outlet is skipped here as
% already set to 1 in above and doesnt matter for other points
inbasin=ismember(order_all,find(fdir<9 & flowdist>=StartDis & flowdist-StartDis <= dowin));
order=order_all(inbasin);

% ID start point at outlet and end point at end of order vector
N = numel(order); %nr*nc;
firstk = find(order==outlet); % position within 'order' where outlet is. 

%StartDis = Dis_order(firstk);

pold=0; %counter

%Only loop through cells w acc<outlet, i.e. upstream
for k = firstk : N
    i = order(k);
    
    %ignore k if it is beyond search radius or downstream of outlet
    %if flowdist(i) < StartDis || flowdist(i)-StartDis > dowin, continue, end;
    
    %if k==firstk, assert(i==outlet); end;
    if basin(i), continue, end; % Skip outlet location
    
    r = rem(i-1, nr)+1;
    c = fix((i-1)/nr)+1;
    % find neighbour
    d = fdir(i);
    rj = r+drow(d);
    cj = c+dcol(d);
    if rj<1,  continue, end;
    if rj>nr, continue, end;
    if cj<1,  continue, end;
    if cj>nc, continue, end;
    % copy catchment membership status from downstream nbr
    basin(r,c) = basin(rj,cj);
    evalcell(r,c)=true;
    
    
    % report - not needed as not all cells used up really
    p = fix(100*k/N);
    if showplot && p>pold
        fprintf('fastfindupstream_Dis Progress: %d%%.\n', p);
        pold=p;
    end;
end
if showplot
    imagescnan(evalcell);axis image
    title(sprintf("# of cells evaluated: %0.0f",sum(evalcell(:))))
end
