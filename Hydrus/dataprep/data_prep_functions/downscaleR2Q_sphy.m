function [Qout, Rout] = downscaleR2Q(Rin, acc, adir, kx, scaleunit, save2path)
% Routes same/low res runoff to high res Q. Res of Q = per res of acc
% Code based on downscaleQ_NAM_loop.m from Gernaat 2017 combined w SPHY
% Currently wc not supported
% Rin - needs to be a a cell array with 2d runoff maps stacked based on timeseries
% acc  - high res flow accumulation map where invalid cells have nan
% adir - high res alternative flow direction map based on Patrick dir conventions
% scaleunit - optional conversion factor to scale runoff depths to volumes 
% kx - flow recession coefficient that accounts for delay in
%       rainfall-runoff response time (0=fast, 1=slow) [see eqn 51 in Terink et al (2015]
% save2path - optional full filepath to save Qout
% Extents on Rin, acc, adir need to be the same
% Qout and Rout are in same units and same extent as acc
% Workflow
%   o	Resize low res to high res using resizer()
%   o	Loop from low to high acc cells to get discharge at cell as sum of local runoff and discharge at all neighbouring cells that drain to it.


% Load directional metrics for Patrick dirs
run(fullfile(pwd,"Hydrus","scripts","defdirs.m" ))

if ~exist('scaleunit')
    scaleunit=1;
end

% Create sorting index to go from  lowest to highest flow accumulation cell
disp('Sorting');
[~, order] = sort(acc(:), 'ascend','MissingPlacement','last'); %Nans sent to last
%nValid=sum(~isnan(acc(:)));

[nr, nc] = size(acc);
N = nr*nc;
nTS=length(Rin);
pold=0;

Qout = cell(1,nTS);%, 'single');  % no wc subtract
Rout = cell(1,nTS);%, 'single');  % no wc subtract

%% Loop through time steps
for mm=1:nTS
    clear Rtmp R Q Qwc
    
    %% Resize R to res of acc if resolution is different
    if all(size(Rin{mm})==size(acc))
        Rtmp = Rin{mm};
    else 
        Rtmp = resizer(Rin{mm},acc,0);
    end
    
    %% convert units
    R=double(Rtmp)*scaleunit;

    % Create placeholder maps
    Q = zeros(nr,nc);%, 'single');  % no wc subtract
    %Qwc = zeros(nr,nc);%, 'single');% wc subtract
    
    %% Compute discharge by downstream routing
    % Loop from lowest to highest acc cells, only valid ones
    for k = 1:N
        % Get address for current cell
        i = order(k);
        r = rem(i-1, nr)+1; % Much faster than [r,c] = ind2sub(sz, i);
        c = fix((i-1)/nr)+1;
        
        % Get local discharge = local runoff
        q = R(r,c);
        
        % Evaluate discharge at cell = local R+ sum of Q from all neighbors that drain to it
        for d = 1:8
            %Get neighbour in d direction
            rj = r+drow(d);
            cj = c+dcol(d);
            % Skip if neighbour is outside domain
            if rj<1,  continue, end
            if rj>nr, continue, end
            if cj<1,  continue, end
            if cj>nc, continue, end
            % Skip neighbour if it does not drain to current cell
            if adir(rj,cj) ~= d, continue, end
            %assert(Q(rj,cj)>0);
            q = q + Q(rj,cj);
        end
        Q(r,c) = q;
        %Qwc(r,c) = q2;
        
        %         If cell is highest acc idx, subtract wcclr as far as possible
        %         if sum(ismember(ida,i))
        %             Qwc(r,c) = max(0,Qwc(r,c) + wcclr(r,c)); %wcclr is negative val or 0
        %         end
    end
    
    % Compiled outputs
    Qout{mm}=Q;
    Rout{mm}=R;
    
    p = fix(100*mm/nTS);
    if p>pold
        fprintf('Timestep %d of %d completed. Progress: %d%%.\n',mm, nTS,p);
        pold=p+0;
    end
end

%% Implement kx shifting similar to SPHY
if ~exist('kx')
    for t=1:12
        if t==1
           Qdelay{t}=(1-kx)*Qout{t}+ kx*Qout{12}; % for jan take leftover flow from dec
        else
        Qdelay{t}=(1-kx)*Qout{t}+ kx*Qout{t-1};
        end
        end
end

%% Save individual mat files for each month
if exist('save2path')
    disp('Saving data')
    matfile = fullfile(save2path, sprintf('Qm%d_wc.mat',mm));
    save(matfile,'-v7.3','Q');
end
