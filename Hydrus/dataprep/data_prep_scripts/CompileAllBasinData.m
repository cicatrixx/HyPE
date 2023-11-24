%% Basin_selector.m
% Prepares the Basin.mat input files with all datasets cropped to basin
% extents that have been pre-processed into .mat files
% Specifically code:
% -Gets extents of each basin in continent based on basin number
% -Extract all cells between identified rectangular extents for 45 active datasets
% -Set all cells falling outside the basin but in rectangular extent to 0/1/9/NaN as per David's conventions
% -Save all basin datsets into one .mat file

% Call from GitConnect
clc
close all
clear all
recheckdata=0;
cp2HPmodel=1;

basin_in = 'UI';  % {'ASIA';'AUS';'CAM';'EUR';'NAM';'SAM';'AFR';'ASIA'};

fprintf('Prepare data for %s \n',basin_in)
root = fullfile(pwd,'data',basin_in,'data'); % get path data folder
minwin=0; % temporary can remove from Hydrus

%forSaving
nbasin=101; HPcurrent='ExistingOnly_3.mat';%5/500/555 David model 1-4 SD model dev 101 SD model final
%nbasin=102; HPcurrent='Existing+UC_3.mat';
nbasin=103; % 103 is model w Qwc that considers irrig consumption

matpath=fullfile(pwd, 'data', 'ASIA', 'Basin_UIB');
matfile = fullfile(matpath,sprintf('PantpeBasin_%d.mat', nbasin));

dnum=0;
%% Load pre-compiled DEM based outputs from theoretical potential run
cellsz_m = 500; % input data cell resolution in meters

%
disp([num2str(dnum),': Reading basin outline']);
fname = fullfile(root,'UI500m_ArcGIS.mat');
basin= load(fname,'catchments','basinmask','acc','dem','fdir','adir', 'flowdist','channel_main_trib', 'channel_ord', 'outside');

if nbasin>=03 %for newer runs select channel definition based on Q
    basin2=     load(fullfile(root,'channel_Qbased.mat'),'channel','channel_main_trib','channel_ord') ;
    basin.channel=basin2.channel;
    basin.channel_main_trib=basin2.channel_main_trib;
    basin.channel_ord=basin2.channel_ord;
end
%
basin_data = basin.basinmask*nbasin;
outside=basin.outside;
%
dnum=dnum+1;disp([num2str(dnum),': Reading subbasin delineation']);
subbasin_data = basin.catchments;
channel_main_trib=basin.channel_main_trib;
channel_ord=basin.channel_ord;
%mask out the misc and endorheic part
basin_data(subbasin_data>=109)=0;
subbasin_data(subbasin_data>=109)=0;
Regions_data=subbasin_data; % this var is used to indicate Image regions, now used as subcatchments
%
dnum=dnum+1;disp([num2str(dnum),': Reading acc']);
acc_data = basin.acc;
%
dnum=dnum+1;disp([num2str(dnum),': Reading Z']);
Z_data = basin.dem;
%
dnum=dnum+1;disp([num2str(dnum),': Reading fdir']);
fdir_data = basin.fdir;
%
dnum=dnum+1;disp([num2str(dnum),': Reading adir']);
adir_data = basin.adir;
%
dnum=dnum+1;disp([num2str(dnum),': Reading flowdist']);
flowdist_data = basin.flowdist;  % distance in number of cells

%% Load Q datasets originally in m3/s: Natural
design_exceedances=[25, 30, 40, 50, 70, 80, 90];

for Qtype="nat"
    dnum=dnum+1;disp([num2str(dnum),': Reading avg annual Q_',Qtype{:}]);
    fname = fullfile(root,sprintf('Qm13_%s.mat',Qtype));
    Q_data = load(fname);
    Q_data = Q_data.Qm13;
    
    
    for qi=1:length(design_exceedances)
        dnum=dnum+1;
        fprintf("%d: Reading Qdesign%d_%s\n",dnum,design_exceedances(qi), Qtype{:});
        fname =fullfile(root,sprintf('Q%d_%s.mat', design_exceedances(qi), Qtype));
        eval(sprintf("Qdesign%d_data = load(fname);", design_exceedances(qi)));
        eval(sprintf("Qdesign%d_Qdesign = Qdesign%d_data.Qdesign;", design_exceedances(qi), design_exceedances(qi)));
        %
        dnum=dnum+1;
        fprintf("%d: Reading Qdesign%d_LF_%s\n",dnum,design_exceedances(qi), Qtype{:});
        eval(sprintf("Qdesign%d_LF_data = Qdesign%d_data.Qdesign_LF;", design_exceedances(qi), design_exceedances(qi)));
        
        % or avg flow through turbine given Qdesign value
        dnum=dnum+1;
        fprintf("%d: Reading Qdesign%d_mean_%s\n",dnum,design_exceedances(qi), Qtype{:});
        eval(sprintf("Qdesign%d_mean_data = Qdesign%d_data.Qdesign_mean;", design_exceedances(qi), design_exceedances(qi)));
    end
end
%
% all(tQdesign90_Qdesign==Qdesign90_Qdesign,'all')
% all(tQdesign90_LF_data==Qdesign90_LF_data,'all')
% all(tQdesign90_mean_data==Qdesign90_mean_data,'all')

%% Load Q datasets: WC
for Qtype="wc"
    dnum=dnum+1;disp([num2str(dnum),': Reading avg annual Q_',Qtype{:}]);
    fname = fullfile(root,sprintf('Qm13_%s.mat',Qtype));
    Qwc_data = load(fname);
    Qwc_data = Qwc_data.Qm13;
    
    for qi=1:length(design_exceedances)
        dnum=dnum+1;
        fprintf("%d: Reading Qdesign%d_%s\n",dnum,design_exceedances(qi), Qtype{:});
        fname =fullfile(root,sprintf('Q%d_%s.mat', design_exceedances(qi), Qtype));
        eval(sprintf("Qdesign%d_data = load(fname);", design_exceedances(qi)));
        eval(sprintf("Qdesign%d_wc_Qdesign = Qdesign%d_data.Qdesign;", design_exceedances(qi), design_exceedances(qi)));
        %
        dnum=dnum+1;
        fprintf("%d: Reading Qdesign%d_LF_%s\n",dnum,design_exceedances(qi), Qtype{:});
        eval(sprintf("Qdesign%d_wc_LF_data = Qdesign%d_data.Qdesign_LF;", design_exceedances(qi), design_exceedances(qi)));
        
        % or avg flow through turbine given Qdesign value
        dnum=dnum+1;
        fprintf("%d: Reading Qdesign%d_mean_%s\n",dnum,design_exceedances(qi), Qtype{:});
        eval(sprintf("Qdesign%d_wc_mean_data = Qdesign%d_data.Qdesign_mean;", design_exceedances(qi), design_exceedances(qi)));
    end
end
% figure
% subplot(1,3,1),imagescnan(Qdesign90_Qdesign-Qdesign90_wc_Qdesign),colorbar
% subplot(1,3,2),imagescnan(Qdesign50_Qdesign-Qdesign50_wc_Qdesign),colorbar
% subplot(1,3,3),imagescnan(Qdesign30_Qdesign-Qdesign30_wc_Qdesign),colorbar

%% Load existing dams datasets in PCS coordinates
dnum=dnum+1;disp([num2str(dnum),': Reading existing dams dataset']);
%fname = fullfile(root,sprintf('ExistingDams2_%d.mat', 3));
fname = fullfile(root,HPcurrent);
load(fname, 'existing_dams','existing_reservoirs');
% use varnames as in Hydrus
r_damst= existing_dams.r_dams ;
c_damst= existing_dams.c_dams ;
GrandIdx= existing_dams.idx_dams ;
Grandlat= existing_dams.yOut;
Grandlon=existing_dams.xOut ;
NoDamsLand = existing_reservoirs ;

%% Load GDPpc ppp data -- Davids data does not have a value for china!
dnum=dnum+1;disp([num2str(dnum),': Reading GDP dataset']);
fname = fullfile(pwd,'data','data_prep','GDPpc','GDPpc2010IsoCode_Indus.csv');
fileID = fopen(fname);
C = textscan(fileID,'%s %s %s %s %s %s %s %s','Delimiter',',','HeaderLines',1);
fclose(fileID);
for i=1:numel(C{3})
    ISOGDP(i,1) = str2num(C{3}{i})'; %ISO number
%    ISOGDP(i,2) = str2num(C{4}{i})'; %GDP pppc2010
    ISOGDP(i,3) = str2num(C{5}{i})'; %GDP current2010 = this is used for cost functions
end

%% Load tech/econ datasets
dnum=dnum+1;disp([num2str(dnum),': Reading existing water bodies from HydroLakes']);
fname = fullfile(root,'WaterBodies.mat');
WaterBodies_data = load(fname);
%
dnum=dnum+1;disp([num2str(dnum),': Reading existing glaciers']);
fname = fullfile(root,'Glaciers.mat');
Glaciers_data = load(fname);
%
dnum=dnum+1;disp([num2str(dnum),': Reading existing glacial lakes']);
fname = fullfile(root,'GlacialLakes.mat');
GlacialLakes_data = load(fname);
AllWaterBodies_data=WaterBodies_data.data>0 | Glaciers_data.data>0 | GlacialLakes_data.data>0 ;
%
dnum=dnum+1;disp([num2str(dnum),': Reading road distance map']);
fname = fullfile(root,'Dis2Road.mat');
DisRoad_data = load(fname);
DisRoad_data = DisRoad_data.data*cellsz_m/1e3; %in km
%
dnum=dnum+1;disp([num2str(dnum),': Reading settlement distance map']);
fname = fullfile(root,'Dis2Settlement.mat');
DisSettlement_data = load(fname);
DisSettlement_data = DisSettlement_data.data*cellsz_m/1e3; %in km
%
dnum=dnum+1;disp([num2str(dnum),': Reading transmission line distance map']);
fname = fullfile(root,'Dis2Transmission.mat');
DisTransmission_data = load(fname);
DisTransmission_data = DisTransmission_data.data*cellsz_m/1e3; %in km
%
dnum=dnum+1;disp([num2str(dnum),': Reading land acq value map']);
fname = fullfile(root,'LandAcq.mat');
LandAcqVal_data = load(fname, 'LandAcqVal');
LandAcqVal_data = LandAcqVal_data.LandAcqVal;
%
dnum=dnum+1;disp([num2str(dnum),': Reading land agri value map']);
fname = fullfile(root,'LandValue_2001-2010_sumLpjml.mat');
LandAgriVal_data = load(fname, 'LandValue_USD_km2');
LandAgriVal_data = LandAgriVal_data.LandValue_USD_km2;
%
dnum=dnum+1;disp([num2str(dnum),': Reading land tree value map']);
fname = fullfile(root,'TreeCost.mat');
LandTreeVal_data = load(fname, 'TreeCost');
LandTreeVal_data = LandTreeVal_data.TreeCost;
%
dnum=dnum+1;disp([num2str(dnum),': Reading population density in inhabitants/km2 per gridcell']);
fname = fullfile(root,'Popd.mat');
Pop_data = load(fname);
Pop_data = Pop_data.data;
%
dnum=dnum+1;disp([num2str(dnum),': Reading countries']);
fname = fullfile(root,'Countries.mat');
Countries_data = load(fname);
Countries_data = Countries_data.data;

%% Load sustainability datasets
dnum=dnum+1;disp([num2str(dnum),': Reading seismic GSHAP map']);
fname = fullfile(root,'GSHAP.mat');
seismicpga_data = load(fname);
seismicpga_data = seismicpga_data.data;
%
dnum=dnum+1;disp([num2str(dnum),': Reading seismic thrust and fault map']);
fname = fullfile(root,'Seismic_Fault_Thrust.mat');
seismicthrust_data = load(fname);
seismicthrust_data = seismicthrust_data.data;
%
dnum=dnum+1;disp([num2str(dnum),': Reading 5 level updated EQ risk map']);
fname = fullfile(root,'geohazards_compiled.mat');
EQ_data = load(fname,'seismic_5level');
EQ_data = EQ_data.seismic_5level;
load(fname,'geohazard_names', 'nhazardlevel','hazard_lnames')
%
dnum=dnum+1;disp([num2str(dnum),': Reading 5 level updated landslide susceptibility map']);
fname = fullfile(root,'geohazards_compiled.mat');
landslide_data = load(fname,'landslide_5level');
landslide_data = landslide_data.landslide_5level;
%
dnum=dnum+1;disp([num2str(dnum),': Reading GLOF path map']);
fname = fullfile(root,'GLOFpath_2deg_Re2.mat');
load(fname,'glofpath_matrix');
glofpath_data = glofpath_matrix>0; % this has GLOF-IDs
%
% dnum=dnum+1;disp([num2str(dnum),': Reading 5 level multi-hazard map']);
% fname = fullfile(root,'geohazards_compiled.mat');
% geohazard_multi_data = load(fname,'geohazard_multi_5levels');
% geohazard_multi_data = geohazard_multi_data.geohazard_multi_5levels;
%
dnum=dnum+1;disp([num2str(dnum),': Reading natural heritage map']);
fname = fullfile(root,'WDPA.mat');
WDPA_data = load(fname);
protected_area_data = WDPA_data.data>0;
%
dnum=dnum+1;disp([num2str(dnum),': Reading cultural heritage map']);
fname = fullfile(root,'CulturalHeritage.mat');
Heritages_data = load(fname);
Heritages_data = Heritages_data.data>0; % this has PA-IDs

% Create buffer around cells
%WDPA_PL20_data=createBuffer(WDPA_PL10_data>0,pa_buffer)-WDPA_PL10_data>0;

%% Set up georefrence
dnum=dnum+1;disp([num2str(dnum),': Reading projection information']);
fname = fullfile(root,'pantpe_georef.mat');
load(fname,'Rw500m','proj500m','Rw100m','proj100m');
Rw = Rw500m;
proj =proj500m;

%% Get UIB bounds for plotting
[x,y]=find(~outside);
k1=boundary(x,y,1);
UIBbounds.x=x(k1);
UIBbounds.y=y(k1);

% figure;imagescnan(outside*0);colormap([1 1 1])
% hold on;
% plot(UIBbounds.y,UIBbounds.x,'k-');

%% Load high-res DEM
dnum=dnum+1;disp([num2str(dnum),': Reading high res DEM']);
load(fullfile(root,'Z100m_fill.mat'))
Z_100m=data;
Rw_100m=Rw100;
disp("Finished loading all factors")

%% Extract relevant data cols for the basin from all datasets; set outside cells to diff vals
[nr,nc]=size(basin_data);
firstrow=1;
lastrow=nr;
firstcol=1;
lastcol=nc;

CurrentBasin = basin_data(firstrow:lastrow, firstcol:lastcol);
channel_main_trib(outside)=0;
channel_ord(outside)=0;
%1
Regions = subbasin_data(firstrow:lastrow, firstcol:lastcol);
Regions(outside) = 0;
%2
acc = acc_data(firstrow:lastrow, firstcol:lastcol);
acc(outside) = 0;
%3
Z = Z_data(firstrow:lastrow, firstcol:lastcol);
Z(outside) = 0;
%4
fdir = fdir_data(firstrow:lastrow, firstcol:lastcol);
fdir(outside) = 9;
%5
adir = adir_data(firstrow:lastrow, firstcol:lastcol);
adir(outside) = 9;
%6
flowdist = flowdist_data(firstrow:lastrow, firstcol:lastcol);
flowdist(outside) = 0;
%27
NoDamsLand(outside) = 0;
%7
Q = Q_data(firstrow:lastrow, firstcol:lastcol);
Q(outside) = 0;
%17
Qwc = Qwc_data(firstrow:lastrow, firstcol:lastcol);
Qwc(outside) = 0;

%Qdesign naturals
for qi=1:length(design_exceedances)
    eval(sprintf("Q%ddesign = Qdesign%d_Qdesign(firstrow:lastrow, firstcol:lastcol);", design_exceedances(qi), design_exceedances(qi)));
    eval(sprintf("Q%ddesign(outside) = 0;",design_exceedances(qi)));
    %
    eval(sprintf("Q%ddesign_LF = Qdesign%d_LF_data(firstrow:lastrow, firstcol:lastcol);", design_exceedances(qi), design_exceedances(qi)));
    eval(sprintf("Q%ddesign_LF(outside) = 0;",design_exceedances(qi)));
    %
    eval(sprintf("Q%ddesign_mean = Qdesign%d_mean_data(firstrow:lastrow, firstcol:lastcol);", design_exceedances(qi), design_exceedances(qi)));
    eval(sprintf("Q%ddesign_mean(outside) = 0;",design_exceedances(qi)));
end

% %Qdesign water consumption
for qi=1:length(design_exceedances)
    eval(sprintf("Q%ddesign_wc = Qdesign%d_wc_Qdesign(firstrow:lastrow, firstcol:lastcol);", design_exceedances(qi), design_exceedances(qi)));
    eval(sprintf("Q%ddesign_wc(outside) = 0;",design_exceedances(qi)));
    %
    eval(sprintf("Q%ddesign_LF_wc = Qdesign%d_wc_LF_data(firstrow:lastrow, firstcol:lastcol);", design_exceedances(qi), design_exceedances(qi)));
    eval(sprintf("Q%ddesign_LF_wc(outside) = 0;",design_exceedances(qi)));
    %
    eval(sprintf("Q%ddesign_mean_wc = Qdesign%d_wc_mean_data(firstrow:lastrow, firstcol:lastcol);", design_exceedances(qi), design_exceedances(qi)));
    eval(sprintf("Q%ddesign_mean_wc(outside) = 0;",design_exceedances(qi)));
end

disp("Finished masking general factors")

%% Mask tech/econ
%29, 30, 31
AllWaterBodies=AllWaterBodies_data(firstrow:lastrow, firstcol:lastcol);
AllWaterBodies(outside) = 0;
%32
DisRoad =DisRoad_data;
DisRoad(outside) = 0;
%33
DisSettlement = DisSettlement_data;
DisSettlement(outside) = 0;
%34
DisTransmission = DisTransmission_data;
DisTransmission(outside) = NaN;
%35
LandAcqVal = LandAcqVal_data(firstrow:lastrow, firstcol:lastcol);
LandAcqVal(outside) = 0;
%36
LandAgriVal = LandAgriVal_data(firstrow:lastrow, firstcol:lastcol);
LandAgriVal(outside) = 0;
%37
LandTreeVal = LandTreeVal_data(firstrow:lastrow, firstcol:lastcol);
LandTreeVal(outside) = 0;
%38
Pop = Pop_data(firstrow:lastrow, firstcol:lastcol);
Pop(outside) = 0;
%39
Countries = Countries_data(firstrow:lastrow, firstcol:lastcol);
Countries(outside) = 0;
disp("Finished masking tech/econ factors")

%% Mask sustainable
%40
seismic_pga = seismicpga_data(firstrow:lastrow, firstcol:lastcol);
seismic_pga(outside) = 0;
%41
seismicthrust = seismicthrust_data(firstrow:lastrow, firstcol:lastcol);
seismicthrust(outside) = 0;
%41
seismic_level = EQ_data(firstrow:lastrow, firstcol:lastcol);
seismic_level(outside) = 0;
%42
landslide_level = landslide_data(firstrow:lastrow, firstcol:lastcol);
landslide_level(outside) = 0;
%43
glofpath = glofpath_data(firstrow:lastrow, firstcol:lastcol);
glofpath(outside) = 0;
%44
Heritages_natural = protected_area_data(firstrow:lastrow, firstcol:lastcol);
Heritages_natural(outside) = 0;
%45
Heritages_cultural = Heritages_data(firstrow:lastrow, firstcol:lastcol);
Heritages_cultural(outside) = 0;

disp("Finished masking sustainability factors")

%% Dummy vars kept so input file is still compatible w original model from David
D=zeros(size(basin_data),'int8');  % temporary placeholder for Depth data

dir=load(fullfile(root,'FDIR500m.mat'),'data');
dir=dir.data;
dir(outside)=0;

seismicpik=zeros(size(basin_data),'int8');
seismicnoaa=zeros(size(basin_data),'int8');
WDPA_PL10 = Heritages_natural;
WDPA_PL20 = zeros(size(basin_data),'int8');
Dis=DisTransmission;
LandValue = LandAgriVal;
Qdesign=Q30design;
Qdesign_mean=Q30design_mean;
Qdesign_LF=Q30design_LF;
Qdesign_wc=Q30design_wc;
Qdesign_mean_wc=Q30design_mean_wc;
Qdesign_LF_wc=Q30design_LF_wc;

%% Save all compiled datasets for basin
if ~isfolder(matpath); mkdir(matpath); end
save(matfile,'acc','Z','fdir','adir','flowdist','dir',...%,'Vlake',
    'Z_100m', 'Rw_100m','cellsz_m',...
    'channel_main_trib', 'channel_ord',...
    'Q','Qdesign','Qdesign_mean','Qdesign_LF',...
    'Qwc','Qdesign_wc','Qdesign_mean_wc','Qdesign_LF_wc',...
    'WDPA_PL10','WDPA_PL20','Heritages_natural','Heritages_cultural',...
    'Regions','Countries','D',...
    'Rw','proj','outside','minwin',...
    'NoDamsLand','r_damst','c_damst','GrandIdx','Grandlat','Grandlon',...
    'AllWaterBodies','LandValue','LandAcqVal','LandAgriVal','LandTreeVal',...
    'Pop','ISOGDP',...
    'Dis','DisRoad','DisSettlement','DisTransmission',...
    'seismicpik','seismicnoaa','seismic_level', 'seismicthrust','landslide_level','glofpath',...
    'geohazard_names', 'nhazardlevel','hazard_lnames',...
    'UIBbounds')
save(matfile,'-regexp',"Q..design", '-append');
%     'Q30design', 'Q30design_mean', 'Q30design_LF',...
%     'Q30design_wc', 'Q30design_mean_wc', 'Q30design_LF_wc',...
%     'Q50design', 'Q50design_mean', 'Q50design_LF',...
%     'Q50design_wc', 'Q50design_mean_wc', 'Q50design_LF_wc',...
%     'Q70design', 'Q70design_mean', 'Q70design_LF',...
%     'Q70design_wc', 'Q70design_mean_wc', 'Q70design_LF_wc',...
%     'Q90design', 'Q90design_mean', 'Q90design_LF',...
%     'Q90design_wc', 'Q90design_mean_wc', 'Q90design_LF_wc',...

%...
disp('All datasets saved')

%% Copy to HPmodel
if cp2HPmodel
    copyfile(matfile,"G:\SurfDrive\HPmodel\data\ASIA\Basin_UIB")
     %  copyfile(matfile,"D:\HPmodel\data\ASIA\Basin_UIB")
    disp("File copied to HPmodel")
end
%% Check if datasets compiled correctly
if recheckdata
    %clc
    clearvars -except matfile
    data=load(matfile);
    namesdata=fieldnames(data);
    ndata=length(namesdata);
    fprintf("No. of vars compiled: %d \n",ndata)
    
    % display arrays as figs and others in command window
    for i=1:ndata
        seld=data.(namesdata{i});
        [r,c] = size(seld);
        if r>1 && c>1
            figure;imagescnan(seld);colorbar;
            title(namesdata{i})
        else
            disp(namesdata{i})
            disp(seld)
        end
    end
end