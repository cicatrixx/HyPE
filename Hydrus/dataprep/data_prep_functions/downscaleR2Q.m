function [Qout, Rout] = downscaleR2Q(Rin, acc, adir, scaleunit, Wcct, verbose, save2path)
% Routes same/low res runoff to high res Q. Res of Q = per res of acc
% Code based on downscaleQ_NAM_loop.m from Gernaat 2017
% Currently wc not supported
% Rin - needs to be a cell array with 2d runoff maps stacked based on timeseries
% acc  - high res flow accumulation map where invalid cells have nan
% adir - high res alternative flow direction map based on Patrick dir conventions
% scaleunit - optional conversion factor to scale runoff depths to volumes
% save2path - optional full filepath to save Qout
% Extents on Rin, acc, adir need to be the same
% Qout and Rout are in same units and same extent as Rin
% Workflow
%   o	Resize low res to high res using resizer()
%   o	Loop from low to high acc cells to get discharge at cell as sum of local runoff and discharge at all neighbouring cells that drain to it.


% Load directional metrics for Patrick dirs
%run(fullfile(pwd,"Hydrus","scripts","defdirs.m" ))
defdirs

if ~exist('scaleunit')
    scaleunit=1;
end

if ~exist('verbose')
    verbose=0;
end

% Create sorting index to go from  lowest to highest flow accumulation cell
disp('Sorting...');
[~, order] = sort(acc(:), 'ascend','MissingPlacement','last'); %Nans sent to last
%nValid=sum(~isnan(acc(:)));

[nr, nc] = size(acc);
N = nr*nc;
nTS=length(Rin);
pold=0;

Qout = cell(1,nTS);%, 'single');  
Rout = cell(1,nTS);%, 'single'); 
Rwcc = cell(1,nTS);%, 'single'); 
Lwcc = cell(1,nTS);%, 'single');  

%% Process water consumption if exists
if exist('Wcct','var')
    % For each cell deduct wc from available runoff
    disp('Subtract water consumption locally')
    wc=1;
    for m=1:nTS
        % Intialize
        Rwcc{m}=zeros(size(Rin{13}),'single'); % Remaining runoff
        Lwcc{m}=zeros(size(Rin{13}),'single'); % Leftover
        
        % Track runoff and leftover after deducting wc
        diffRW=Rin{m}-Wcct{m};
        Rwcc{m}(diffRW>=0)=diffRW(diffRW>=0); % remaining runoff after local water consumption
        Lwcc{m}(diffRW<0)=diffRW(diffRW<0); % left over water demand that need to be met
    end
    
    % Find highest acc idx in 0.5x0.5 degree blocks
    % For cells w unmet demands, need to subtract remaining water consumption from highest flow accumulation cell in routing process
    [~,~,ida] = MaxIDBlockFinder(acc);
    sub_wcl_cells= acc*0;
    sub_wcl_cells(ida)=1;
    
else
    Rwcc=Rin;
    wc=0;
    %     for m=1:nTS
    %         Wcct{m}=0*Rin{m};
    %     end
end


%% Loop through time steps
disp('Routing discharge...')

for mm=1:nTS
    clear R Q
    %% Resize R to res of acc and convert units from mm to m3
    R=double(resizer(Rwcc{mm},acc,0))*scaleunit; %need to do this conversion to double to not lose precision
    if wc==1
        wcl=double(resizer(Lwcc{mm},acc,0))*scaleunit; %double(resizer(Wct{mm},acc,0))*scaleunit; %need to do this conversion to double to not lose precision
    end
    
    % Create placeholder maps
    Q = zeros(nr,nc);%, 'single');
    
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
        if wc==1 && sub_wcl_cells(r,c)==1
            Q(r,c) = max(0,Q(r,c) - wcl(r,c)); 
        end
    end
    
    % Compiled outputs
    Qout{mm}=Q;
    Rout{mm}=R;
    
    p = fix(100*mm/nTS);
    if verbose && p>pold
        fprintf('Timestep %d of %d completed. Progress: %d%%.\n',mm, nTS,p);
        pold=p+0;
    end
end

%% Save each time step separately
if exist('save2path')
    disp('Saving data')
    matfile = fullfile(save2path, sprintf('Qm%d_wc.mat',mm));
    save(matfile,'-v7.3','Q');
end
