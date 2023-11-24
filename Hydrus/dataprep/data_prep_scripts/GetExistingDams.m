% Load all dams and select only ones in basin. Project data if not in PCS.
% Shift point to nearest neighbour with Qneighbour>2*Qfocus. Delineate
% reservoir or NoDamLand using FFupstream. If dam height is available, use it to cutoff
% upstream else just use cell location.
% Save only existing and existing+UC dams for use in model. did not filter
% by storage or other plant type cos if its there right now, it takes up
% the land for new projects regardless!

% For GCS, xOut=lonOut; yOut=latOut
% For PCS, additonally need a proj input

% All grandlocs are assumed to be dam locations
% Z used to get elevation of outlets, (Zoutlets(k)+grand_height(k)) used to
% get reservoir above dam outlet where height is available
% Flowdist, acc and fdir used for tracing upstream using FFupstream
% Q used to get RCs and shift dam outlets to
% higher Q

% Merged david's findGrandDams and NoDamsLand functions into one to
% pre-process Grand dams related inputs
% Found 11 dams in UI in GCS and 12 in PCS right at the border

clear all; %close all; clc

root=pwd; %%G:/SurfDrive/GitConnect/
showcheckswDavid=0; % compare w Davids outputs
showplots=0;

%res_sradius = 10000/500; %~assumes reservoir cannot be more than 50km upstream

% Input datasets
nbasin=3; %res=500; % in m
basindata_fullpath=fullfile(root,sprintf('data/ASIA/Basin_UIB/PantpeBasin_%d.mat', nbasin));
load(basindata_fullpath,'Z','Q','flowdist','acc','fdir','adir','Rw','outside')
coords=Rw.CoordinateSystemType;

%SD visualized DB - Load all plants w latlon info
visualizeddb_fullpath=fullfile(pwd,'data/data_prep/Grand/Clean_Mix_GIS_allHPwLatLon.xlsx'); %%
dbcoords="geographic";

%GrandDB from David
%existingdb_fullpath=fullfile('G:/SurfDrive/GitConnect/data/data_prep/Grand/grand_dams_latlonDH.csv');
% dbcoords="geographic"
% existingdb_fullpath=fullfile('G:/SurfDrive/GitConnect/data/data_prep/Grand/grand_dams_PCS_DH.csv');
% dbcoords="planar";

% get defdirs
defdirs
%load('G:\SurfDrive\GitConnect\data\ASIA\Basin\defdirs.mat','drow','dcol')

% Outlocations
matpath=fullfile(pwd,'data/UI/data');
outtype="Existing+UC+All"; %"Existing+UC"; %"ExistingOnly"; %  
outmat_fullpath=fullfile(matpath,sprintf('%s_%d.mat', outtype, nbasin));
outtif_allres  =fullfile(matpath,sprintf('NoDamsLand_ALLReservoirs_IDs_%d.tif',nbasin));
outtif_allres_stats  =fullfile(matpath,sprintf('NoDamsLand_ALL_Reservoirs_Status_%d.tif',nbasin));
%outtif_currentres=fullfile(matpath,sprintf('NoDamsLand_CurrentReservoirs_%d.tif',nbasin));
fprintf("\n\nLoaded inputs for Basin_%d in coords %s\n", nbasin,coords)

%% Read in existingDB for the globe
% Missing data is -99
dams_all = readtable(visualizeddb_fullpath); %
dams=dams_all{:,1:6};

%% Convert locations to r,c for dams within the basin latlon limits
if strcmp(coords,'geographic') && strcmp(dbcoords,'geographic')
    seldams = dams(:,5) > Rw.Lonlim(1) & dams(:,5) < Rw.Lonlim(2) & dams(:,6) > Rw.Latlim(1) & dams(:,6) < Rw.Latlim(2);
    [r_dams1, c_dams1] = setpostn(Q, Rw, dams(seldams,6), dams(seldams,5)); %Convert latitude-longitude to data grid rows and columns
elseif strcmp(coords,'planar') && strcmp(dbcoords,'geographic')
    % First project to PCS and then crop
    load(basindata_fullpath,'proj')
    [dams_X, dams_Y] = projfwd(proj,dams(:,6), dams(:,5));
    seldams = dams_X> Rw.XWorldLimits(1) & dams_X < Rw.XWorldLimits(2) & dams_Y > Rw.YWorldLimits(1) & dams_Y < Rw.YWorldLimits(2);
    [r_dams1, c_dams1] = worldToDiscrete(Rw, dams_X(seldams), dams_Y(seldams)); %Convert XY to data grid rows and columns
    % confirmed that this matched Arcgis output
elseif strcmp(coords,'planar') && strcmp(dbcoords,'planar')
    dams_X = dams(:,5);
    dams_Y = dams(:,6);
    seldams = dams_X> Rw.XWorldLimits(1) & dams_X < Rw.XWorldLimits(2) & dams_Y > Rw.YWorldLimits(1) & dams_Y < Rw.YWorldLimits(2);
    [r_dams1, c_dams1] = worldToDiscrete(Rw, dams_X(seldams), dams_Y(seldams)); %Convert XY to data grid rows and columns
    % yet to confirm that this matches rasterization from gdal
end
db_id=dams(seldams,1);

fprintf("No. of dams within bounds: %d\n",numel(db_id))

%% Get  dam details if exists in current basin
if isempty(r_dams1)
    %% If no valid VisualizedDams in current basin, set outputs to zero and close function
    NoDamLand = zeros(size(Q),'int8');
    db_id=-99;r_dams_shift=-99;c_dams_shift=-99;idx_dams_shift=-99;latOut=-99;lonOut=-99;xOut_shifted=-99;yOut_shifted=-99;
    shiftedrc=-99;dam_height=-99; lake_cellcnt=-99;
    disp('No dams in current basin')
else
    %% If valid VisualizedDams, select only dams inside the basin domain (previously this was based on Q=0)
    idx_dams = sub2ind(size(Z),r_dams1,c_dams1); % dam index in full data window
    selitem=outside(idx_dams)~=1; %sel dams indexes not =1 in outside matrix
    idx_dams=idx_dams(selitem);
    r_dams=r_dams1(selitem);
    c_dams=c_dams1(selitem);
    db_id=db_id(selitem);
    fprintf("No. of dams within basin: %d\n",numel(db_id))
    
    % Get dam details
      % Get elevation at each dam bottom and top for dam_ids that match select dams
    [~,LocInDams]= ismember(db_id,dams(:,1));
    dam_height = dams(LocInDams,2); %Heights
    energy_GWh=dams(LocInDams,3);
    energy_MW=dams(LocInDams,4);
    status=dams_all{LocInDams,7};
    planttype=dams_all{LocInDams,9};
    % Throw error if any dam has Q=0; should probably set variable Qshift_win for these manually
    assert(~any(Q(idx_dams)==0))
    
    %% Shift VisualizedDB locations to nearest neighbours cells with Qneighbour>2*Qfocus
    disp('Shifting VisualizedDB locations to higher Q points')
    Qshift_win=1;
    % Initialize
    Qdams_preCorr = Q(sub2ind(size(Q),r_dams,c_dams));
    r_dams_shift = r_dams;
    c_dams_shift = c_dams;
    shiftedrc=logical(0*c_dams_shift); % Tracker for dams that were shifted
    
    for i=1:numel(r_dams)
        % Create a window of one cell around the current cell, i.e. nearest neighbor
        firstrow = r_dams(i)-Qshift_win;
        lastrow  = r_dams(i)+Qshift_win;
        firstcol = c_dams(i)-Qshift_win;
        lastcol  = c_dams(i)+Qshift_win;
        Qwin = Q(firstrow:lastrow,firstcol:lastcol); %Select small Q window
        
        Qfocus = Q(r_dams(i),c_dams(i));
        Qmax = max(Qwin(:));
        Qratio = Qmax/Qfocus;
        
        if Qratio>2
            shiftedrc(i)=1;
            [r_max,c_max]=find(Qwin==Qmax); %Find max
            
            %Shift factor to move focus r,c to neighbouring highest point in window grid
            r_trans = r_max(1)-Qshift_win-1; %In case multiple max, pick first option
            c_trans= c_max(1)-Qshift_win-1;
            
            %Transform to large grid
            r_dams_shift(i) = r_dams(i)+r_trans;
            c_dams_shift(i) = c_dams(i)+c_trans;
            
            %check
            if showplots
                figure;clf;imagesc(Qwin);axis image;colormap(flipud(gray)); hold on
                plot(Qshift_win+1,Qshift_win+1,'r.'); % Original focus point
                plot(c_max(1),r_max(1),'*b','markersize',20)        % Shifted focus point
                %hold off
            end
        end
    end
    % Update indices
    idx_dams_shift=sub2ind(size(Z),r_dams_shift,c_dams_shift);
    fprintf('No. of dams shifted: %d\n',sum(shiftedrc));
    
    %% Get new coordinates after shifting
    if strcmp(coords,'geographic')
        [latOut, lonOut] = setltln(Q, Rw, r_dams_shift, c_dams_shift); %Coordinates Visualized dams
        xOut_shifted=lonOut; yOut_shifted=latOut;
        
    elseif strcmp(coords,'planar')
        [xOut_shifted, yOut_shifted] = intrinsicToWorld(Rw, c_dams_shift, r_dams_shift) ;
    end
    
    %% Double check w David's for 15sUI GCS basin500
    if showcheckswDavid && strcmp(coords,'geographic')
        disp("Getting david's grand dam shift")
        
        orig=load(fullfile(root,'data/ASIA/Basin/Basin_500.mat'),'Q','Rw');
        root_bil = 'G:\SurfDrive\GitConnect';
        [GrandR,GrandC,GrandIdx,Grandlat,Grandlon] = findGrandDams(root_bil,orig.Rw,orig.Q);
        
        all(r_dams_shift==GrandR' & c_dams_shift==GrandC'& idx_dams_shift==GrandIdx' & lonOut==Grandlon' & latOut==Grandlat')
        
        figure;clf;imagesc(Q); axis image;colormap(flipud(gray));hold all
        set(gca, 'ColorScale', 'log')
        hold on;
        plot(c_dams1,r_dams1,'.b','markersize',15,'DisplayName','Original in Bounds'); %Existing dams
        legend('-DynamicLegend')
        plot(c_dams,r_dams,'ok','markersize',20,'DisplayName','Original in Basin');
        plot(c_dams_shift(shiftedrc),r_dams_shift(shiftedrc),'+r','markersize',10,'DisplayName','Shifted');
        plot(GrandC,GrandR,'oc','markersize',15,'DisplayName','Davids Vals'); %Existing dams
    end
    
    %% Get NoDamsLand - reservoirs upstream of dams
    %Initialize
    lake_cellcnt = zeros(size(r_dams_shift));
    all_reservoirs = nan(size(Q));
        
    Zdambottom = single(Z(idx_dams_shift));
    Zdamtop=single(Zdambottom+dam_height);
    fprintf("No. of dams w no dam height info: %d of %d\n", sum(dam_height==-99),numel(dam_height))
    
    % if dam height available, cut off upstream to dam elevation, else take
    % upstream upto search radius
    disp('Delineating dam reservoirs upstream')
    
    for k=1:numel(r_dams_shift)
        %fprintf('Outlet #%d of %d\n',k,numel(r_dams_shift))
        clear RDupstream Damlake
        
        if dam_height(k)==-99 %If DB gives nodata, consider plant has no reservoir and only use the cell location as covered
            Damlake=logical(size(Q)); %fastfindupstream_DisSD4(acc,fdir,adir,flowdist,drow,dcol,idx_dams_shift(k),0, res_sradius);
            Damlake(idx_dams(k)) =1;
        else
            RDupstream = fastfindupstream_DisSD4(acc,fdir,adir,flowdist,drow,dcol,idx_dams_shift(k),0);
            Damlake = RDupstream & Z < Zdamtop(k);
        end
        
        lake_cellcnt(k) = sum(Damlake,'all');
        all_reservoirs(Damlake)=db_id(k);        
    end
    %% Reclassify reservoirs based on status % 50=Existing, 20=Under-Construction, 1 =Raw or planned
    all_reservoirs_status=nan(size(all_reservoirs));
for k=1:height(status)
    if strcmp(status(k),"Existing")
        all_reservoirs_status(all_reservoirs==db_id(k)) = 50;
    elseif strcmp(status(k),"Under Construction")
        all_reservoirs_status(all_reservoirs==db_id(k)) = 20;
    else
        all_reservoirs_status(all_reservoirs==db_id(k)) = 1;
    end
end

    %% Compare w original
    if showcheckswDavid
        if strcmp(coords,'geographic')
            disp("Getting david's no dam land")
            orig=load(fullfile(root,'data/ASIA/Basin/Basin_500.mat'),'Z','Q','flowdist','Dis','D','acc','fdir','minwin');
            load(fullfile(root,'data/ASIA/basins.mat'),'basin');
            orig.basin=basin;
            origNoDamLand=NoDamsLand(pwd,pwd,'ASIA',5,orig.Z,orig.Q,orig.flowdist,orig.Dis,orig.D,orig.acc,single(orig.fdir),orig.minwin,orig.basin);
            all(origNoDamLand(:)==all_reservoirs(:))
            figure;clf;imagesc(origNoDamLand==all_reservoirs); axis image
            title("Compare orig=current nodamsland")
            
        elseif strcmp(coords,'planar')
            load('G:\SurfDrive\GitConnect\data\UI\data\WaterBodies.mat', 'data')
            figure;clf;imagesc(single(data>0)+single(all_reservoirs)*10); axis image;%colormap(flipud(gray));hold all
            title(" nodamsland + HydroLAKES")
        end
    end
        
    %% Plot sites before and after w nodamsland
    figure;clf;imagescnan(maskBasin(all_reservoirs,~outside)); axis image;colormap(flipud(gray));hold all
    hold on;
    % plot(c_dams1,r_dams1,'.b','markersize',15,'DisplayName','Original in Bounds'); %Existing dams
    legend('-DynamicLegend')
    plot(c_dams,r_dams,'.b','markersize',5,'DisplayName','Original in Basin');
    plot(c_dams_shift(shiftedrc),r_dams_shift(shiftedrc),'+r','markersize',10,'DisplayName','Shifted');
    title(sprintf("Initial and final locations of ALL dams and reservoirs in Basin_%d in coords %s\n", nbasin,coords))
    
    %% Create Q overlaid w new reservoirs
    tmp=all_reservoirs;
    tmp(isnan(tmp))=0;
    Q_img= truecolorsc(log(Q),flipud(gray));
    Q_img = burnmask(Q_img, ~tmp);
    figure
    imagesc(Q_img); axis image;
    hold on;
    % plot(c_dams1,r_dams1,'.b','markersize',15,'DisplayName','Original in Bounds'); %Existing dams
    legend('-DynamicLegend')
    plot(c_dams,r_dams,'.b','markersize',5,'DisplayName','Original');
    plot(c_dams_shift(shiftedrc),r_dams_shift(shiftedrc),'+r','markersize',10,'DisplayName','Shifted');
    text(c_dams,r_dams,categorical(db_id),'FontSize' ,6,'Color','g')%,'BackgroundColor',.9*[1 1 1]);%,'FontWeight','bold'
    
    %%  Select only current(existing or under construction) dams
    if outtype=="Existing+UC"
        current_rows= (dams_all{:,7}=="Existing") | (dams_all{:,7}=="Under Construction");
    elseif outtype=="ExistingOnly"
        current_rows= (dams_all{:,7}=="Existing") ;
    elseif outtype=="Existing+UC+All"
        current_rows=1:numel(db_id);
    end
    current_id = db_id(current_rows);
    current_res =zeros(size(Q));
    for k=current_id'
        current_res(all_reservoirs==k)= k;
    end
    
    Q_img= truecolorsc(log(Q),flipud(gray));
    
    figure;clf;
    imagescnan(maskBasin(current_res,~outside));
    imagescnan(burnmask(Q_img, ~current_res));
    %colormap(flipud(gray))
    axis image;hold all
    hold on;
    % plot(c_dams1,r_dams1,'.b','markersize',15,'DisplayName','Original in Bounds'); %Existing dams
    legend('-DynamicLegend')
    plot(c_dams_shift(current_rows),r_dams_shift(current_rows),'+r','markersize',9,'DisplayName','Current HPs');
    scatter(c_dams_shift(current_rows),r_dams_shift(current_rows),rescale(energy_GWh(current_rows),1,1000),'DisplayName','Current HPs')
    title(sprintf("Current hydropower plant and their tentative reservoirs in coords %s",coords))     
    end

%% Save to matfile
all_HPs = table(db_id, energy_GWh, energy_MW, status, r_dams_shift, c_dams_shift, idx_dams_shift, xOut_shifted, yOut_shifted, shiftedrc,dam_height,lake_cellcnt,planttype,...
    'VariableNames',{'id_dams', 'GWh', 'MW','status', 'r_dams','c_dams', 'idx_dams','xOut','yOut','shiftedrc','h_dams','cellcnt_dams','planttype'});
existing_dams = all_HPs(current_rows,:);
existing_reservoirs = current_res;
if exist ('outmat_fullpath')
    save(outmat_fullpath,'existing_dams', 'existing_reservoirs','nbasin','coords','all_reservoirs_status')
    fprintf('Saved output to %s\n',outmat_fullpath)
    % Save visualized data to excel
    scfname=strrep(outmat_fullpath,".mat",".xlsx");
    writetable(existing_dams,scfname,'Sheet',"Visualized")
end

%% Save nodamsland to tif
setnodatavals=-9999;
if exist('outtif_allres') && strcmp(coords,'geographic')
    geotiffwrite(fullfile(grandloc,outtif),all_reservoirs,Rw)
elseif exist('outtif_allres') && strcmp(coords,'planar')
    all_reservoirs(isnan(all_reservoirs))=setnodatavals;
        all_reservoirs_status(isnan(all_reservoirs_status))=setnodatavals;
geotiffwrite(outtif_allres,all_reservoirs,Rw,'GeoKeyDirectoryTag',proj.GeoTIFFTags.GeoKeyDirectoryTag)
    geotiffwrite(outtif_allres_stats,all_reservoirs_status,Rw,'GeoKeyDirectoryTag',proj.GeoTIFFTags.GeoKeyDirectoryTag)
        %geotiffwrite(outtif_currentres,current_res,Rw,'GeoKeyDirectoryTag',proj.GeoTIFFTags.GeoKeyDirectoryTag)
    fprintf('Saved output to %s\n',outtif_allres)
end