function [lrwidth_out genwidth_out]=AsymDamLenght(root,continent_in,lat,lon,winsize,RDZdiff,ASYMDP)

% clear all
% root = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global';
% continent_in = 'EUR';
% lat=lat{k}{l}(m);
% lon=lon{k}{l}(m);
% winsize=30;
% RDZdiff = 100;
% ASYMDP=1;

%% Aproximate grid cell size
dx = 90 * cosd(lat);   % shrinks towards poles
dy = 90;               % constant
dxy = sqrt(dx^2+dy^2); % diagonal

%% Load hi-res (3'') Hydrosheds DEM and DIR

dem = hs3s_window('con',lat,lon, winsize, root,continent_in);

if dem==0;
    lrwidth_out=NaN;
    genwidth_out=NaN;
    return;
end;

% Convert from int->single and set NaN values for ocean cells
ocean = dem==-9999;
dem = single(dem);
dem(ocean) = NaN;

ldd = hs3s_window('dir',lat,lon, winsize, root,continent_in);
lddvals = unique(ldd(:))';

if ldd==0;
    
    disp('No dir3s')
    lrwidth_out=NaN;
    genwidth_out=NaN;
    
%     [nrows,ncols] = size(dem);
%     DX = dx * abs((1:ncols) - winsize - 1);
%     DY = dy * abs((1:nrows) - winsize - 1);
%     DX = repmat(DX, [nrows 1]);
%     DY = repmat(DY', [1 ncols]);
%     flydist = sqrt(DX.^2 + DY.^2);
%     
%     zsteps = 1:RDZdiff; % elevation steps to consider
%     
%     genwidth = zeros(size(zsteps)); % corresponding generic width
%     
%     for i=1:numel(zsteps);
%         % generic width
%         idx = dem>zsteps(i);
%         if isempty(min(flydist(idx)))==1;
%             genwidth(i)=max(genwidth);
%         else
%             genwidth(i) = 2*min(flydist(idx));
%         end
%     end
%     
%     genwidth(find(~genwidth))=90;
%     
%     genwidth_out=genwidth(end);
%     lrwidth_out=genwidth_out;
    return;
end;

% Convert from ESRI formatt, i.e. 1,2,4,8,...128 for E, SE, S etc. to
% consecutive numbers 1,2,3,..8.

defdirs_asym;

% Convert
ldd2 = ldd;
ldd2(ldd==0)   = OUTLET;
ldd2(ldd==1)   = EAST;
ldd2(ldd==2)   = SEAST;
ldd2(ldd==4)   = SOUTH;
ldd2(ldd==8)   = SWEST;
ldd2(ldd==16)  = WEST;
ldd2(ldd==32)  = NWEST;
ldd2(ldd==64)  = NORTH;
ldd2(ldd==128) = NEAST;
ldd2(ldd==247) = SINK;

% Replace
ldd = ldd2;
clear ldd2;

lddvals = unique(ldd(:))';

%% Compute neighbourhood grid
nbr = nbrtable(dem, 8);

%% Compute map of downstream nbrs
down = zeros(size(ldd), 'int32'); % Create empty map
for i=1:numel(ldd),
    if ldd(i)==OUTLET, continue; end; % Skip outlet locations
    if ldd(i)==SINK, continue; end; % Skip sink points, oceans etc.
    j = nbr(i, ldd(i));             % Get neighbour address
    if j==0, continue, end;         % Neighbour does not exist
    down(i) = j;                    % Assign neighbourhood index
end;

%% Find outlet

rdam = winsize+1; % by definition
cdam = winsize+1;
idam = sub2ind(size(dem), rdam, cdam);

i = idam;
while down(i)>0
    i = down(i);
end;
outlet = i;
[routlet, coutlet] = ind2sub(size(dem), outlet);

% xoutlet = x(coutlet);
% youtlet = y(routlet);

%% Show
% figure(1);
% imagesc(x,y,dem);
% axis equal; axis tight;
% title('Elevation');
% xlabel('Longitude');
% ylabel('Latitude');
%
% % Mark target
% hold on;
% hdam = plot(lon, lat,  'r.', 'markersize',24);
% hold off;
% legend(hdam, 'Dam site', 'Location','BestOutside');
%
% % Mark outlet
% hold on;
% houtlet = plot(xoutlet, youtlet, '.', 'markersize',24);
% set(houtlet, 'color', 0.8*[1 1 1]);
% hold off;
% legend([hdam houtlet], 'Dam','Outlet', 'Location','BestOutside');

%% Contributing area
A = down2carea(nbr,down,dem);

% imagesc(x,y,log(A));
% axis equal; axis tight;
% colormap(flipud(gray(256)));
%
% hold on;
% hdam = plot(lon, lat,  'r.', 'markersize',24);
% houtlet = plot(xoutlet, youtlet, '.', 'markersize',24);
% set(houtlet, 'color', 'b');
% hold off;
%
% xlabel('Longitude');
% ylabel('Latitude');
% title('Contributing area');
% legend([hdam houtlet], 'Dam','Outlet', 'Location','BestOutside');

%% Flow distance map
% first set up a matrix of local flow distances
ldist = zeros(size(ldd));
ldist(ldd==EAST)  = dx;
ldist(ldd==SOUTH) = dy;
ldist(ldd==WEST)  = dx;
ldist(ldd==NORTH) = dy;
ldist(ldd==SEAST) = dxy;
ldist(ldd==SWEST) = dxy;
ldist(ldd==NWEST) = dxy;
ldist(ldd==NEAST) = dxy;

flowdist = zeros(size(dem));
[~,upstream] = sort(A(:), 'descend'); % Move from the outlet upstream
for k=1:numel(dem);
    i = upstream(k);
    if down(i)==0
        if i==outlet
            flowdist(i) = 0;
        else
            flowdist(i) = NaN;
        end
    else
        flowdist(i) = flowdist(down(i))+ldist(i);
    end
end
%
% imagesc(x,y,flowdist);
% axis equal; axis tight;
% colormap default;
%
% hold on;
% hdam = plot(lon, lat,  'r.', 'markersize',24);
% houtlet = plot(xoutlet, youtlet, '.', 'markersize',24);
% set(houtlet, 'color', 0.8*[1 1 1]);
% hold off;
%
% xlabel('Longitude');
% ylabel('Latitude');
% title('Flow distance');
% legend([hdam houtlet], 'Dam','Outlet', 'Location','BestOutside');

%% Search for inlet
% First define the border
border = zeros(size(dem));
border(1,:)   = 1;
border(end,:) = 1;
border(:,1)   = 1;
border(:,end) = 1;

% Then move upstream from the outlet
i = outlet;
going = true;
for k=1:numel(dem)
    count=0;
    for d=1:8
        j = nbr(i,d);
        if j==0, continue; end; % skip non-existing nbrs
        if down(j)~=i, continue; end; % skip non-contributing nbrs
        
        if count==0
            % first serious candidate
            propose = true;
        elseif dem(j) < zcandidate,
            % replace by lower upstream nbr
            propose = true;
        elseif dem(j) == zcandidate && A(j) > Acandidate
            % replace by equal elevation, but larger area
            propose = true;
        else
            propose = false;
        end;
        
        if propose
            dcandidate = d;
            jcandidate = j;
            zcandidate = dem(j);
            Acandidate = A(j);
            count = count+1;
        end;
    end;
    
    i=jcandidate;
    if flowdist(i)>winsize
        % only consider terminating if we have travelled half the grid
        if border(i), break; end; % stop if we hit the border
        if count==0,  break; end; % terminate at source
    end
end

inlet = i;
[rinlet,cinlet] = ind2sub(size(dem), inlet);
% xinlet = x(cinlet);
% yinlet = y(rinlet);

%% Find main stream
% Finding the main stream now is now easy: just trace downstream from the inlet.
% channel cell type constants
CHANNELHEAD = 1;
CHANNEL     = 2;
LINKHEAD    = 3;
LINKTAIL    = 4;

channel = false(size(dem));
ctype   = zeros(size(dem),'uint8');
channel(inlet) = true;
ctype(inlet)   = CHANNELHEAD;
i = inlet;
while down(i)
    i = down(i);
    channel(i) = true;
    ctype(i) =  CHANNEL;
end

% imagesc(x,y,channel);
% axis equal; axis tight;
% xlabel('Longitude');
% ylabel('Latitude');
% title('Main channel');

%% Define the banks
% Bank types
HILLSLOPE  = 1;
CHANNEL    = 2;
SOURCEBANK = 3;
LEFTBANK   = 4;
RIGHTBANK  = 5;
STRANGEBANK= 6;
ERROR      = 7;

bank  = findbanks(A,channel,ldd,nbr,down,ctype);

% create color palette for banks
bankpal = [1 1 1;... % white for non-land areas
    0.5 0.5 0.5;... % gray for non-assigned hillslope cells
    0 0 1;... % blue for channel cells
    1 0 0;... % red for source cells
    0 1 0;... % green for left banks
    1 1 0;... % yellow for right banks
    1 0 1;... % purple for strange banks
    0 0 0];   % black for errors

% image(x,y,bank);
% colormap(bankpal);
% h=colorbar;
% set(h,'ytick', 1.5 : 7.5, ...
%     'yticklabel',{'unassigned','channel','source','left','right','strange','error'});
% axis equal; axis tight;
% xlabel('Longitude');
% ylabel('Latitude');
% title('Bank types');

%% Extend banks to hillslopes
% Next step is to assign the bank type (left,right) to the corresponding
% hillslope (and non-mainstream) cells. If we work our way uphill, we can
% just simply copy the type from the downhill nbr cell.

hstype = bank;
[~,uphill] = sort(A(:), 'descend'); % walk uphill
for k=1:numel(A)
    i = uphill(k);
    if bank(i)==HILLSLOPE
        j = down(i);
        if j>0
            if hstype(j)==SOURCEBANK
                hstype(i) = SOURCEBANK;
            end;
            if hstype(j)==LEFTBANK
                hstype(i) = LEFTBANK;
            end;
            if hstype(j)==RIGHTBANK
                hstype(i) = RIGHTBANK;
            end;
        end
    end
end

% imagesc(hstype);
% colormap(bankpal);
% h=colorbar;
% set(h,'ytick', 1.5 : 7.5, ...
%     'yticklabel',{'unassigned','channel','source','left','right','strange','error'});
% axis equal; axis tight;
% xlabel('Longitude');
% ylabel('Latitude');
% title('Bank types');

%% Distance to focal point
% Compute the distance (in m) to the center point

[nrows,ncols] = size(dem);
DX = dx * abs((1:ncols) - winsize - 1);
DY = dy * abs((1:nrows) - winsize - 1);
DX = repmat(DX, [nrows 1]);
DY = repmat(DY', [1 ncols]);
flydist = sqrt(DX.^2 + DY.^2);

%% Dam height vs width
% create indices of where the banks are
leftbank = hstype==LEFTBANK;
rightbank = hstype==RIGHTBANK;

r0 = winsize+1; % row,column of focal point
c0 = winsize+1;
zshift = dem - dem(r0,c0); % shifted DEM, focal point at z=0

zsteps = 1:RDZdiff; % elevation steps to consider

genwidth = zeros(size(zsteps)); % corresponding generic width
lwidth = zeros(size(zsteps)); % halfwidth (left bank)
rwidth = zeros(size(zsteps)); % halfwidth (right bank)
lrwidth = zeros(size(zsteps)); % combined width (left+right)
%%
for i=1:numel(zsteps);
    % generic width
    idx = dem>zsteps(i);
    if isempty(min(flydist(idx)))==1;
        genwidth(i)=max(genwidth);
    else
        genwidth(i) = 2*min(flydist(idx));
    end
    
    % left halfwidth
    lidx = idx & leftbank;
    if isempty(min(flydist(lidx)))==1;
        lwidth(i)=max(lwidth);
    else
        lwidth(i) = min(flydist(lidx));
    end
    
    % right halfwidth
    ridx = idx & rightbank;
    if isempty(min(flydist(ridx)))==1;
        rwidth(i)=max(rwidth);
    else
        rwidth(i) = min(flydist(ridx));
    end
    
    % left+right width
    lrwidth(i) = lwidth(i) + rwidth(i);
end;

%% Insert minimum dam width of 90m, size of a 3s cell

lrwidth(find(~lrwidth))=90;
genwidth(find(~genwidth))=90;

%%
if ASYMDP==1
    lrwidth_out=lrwidth(end);
    genwidth_out=genwidth(end);
else
    lrwidth_out=lrwidth;
    genwidth_out=genwidth;
end

% clf;
% h=plot(genwidth,zsteps,'k-', lrwidth,zsteps,'r', lwidth,zsteps,'g', rwidth,zsteps,'b');
% set(h(1:2), 'linewidth',1);
% legend('generic','left+right', 'left (half)','right(half)');
% xlabel('Dam width (m)');
% ylabel('Dam height (m)');

end