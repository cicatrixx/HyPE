%% Compile SPHY and PCRASTER based .map files for Q/R and .tiff files for acc, ldd and dem into one .mat file
% Read 500m .map SPHY outputfiles for TotR (for 12 months) and QAll (for 12 months+yr) - 40 yr avg climatology from (1979-2018)
% Generate annual TotR and QAll as M13 by average monthly files. Loads .maps as .tiffs then save into .mat files
% Load .tiff for acc, ldd, dem and basin mask prepared elsewhere

% Filled cells outside basins same way as done by David in Basin_selector
% for consistency as some functions dont work well, esp fdir and R/Q.
% Things are kept in original extents for now. No masking

% USES GDAL_TRANSLATE!
% Created By    : Sanita Dhaubanjar on 28 May 2020
% Created For	: SustaIndus project WP2
%=========================
clc
close all
clear all
%% Inputs
% tiff nan defaults
nodata32bitFlo = -340282346638528859811704183484516925440.00; % for 32 bit signed FLOAT type
nodata32bitInt = -2147483648;
nodata8bit = 255; % for 8 bit unsigned
%
showplot=0;
prepAnnualQ_T=0;
loadR_Q=0;
%output locs
matpath='G:/GitConnect/data/UI/';
%matpath=
QRout='SPHY5km_40yrClimatology.mat';
data500m='UI500m.mat'; % for Z, acc, fdir and adir

%% Filenames - 500m
root500m ='G:/SPHY/Arthur/500m';%pwd;
scn="corrected";
if scn=="corrected"
    demfile = 'dem_500m_v3SD_matlabGDAL.tif';
    dirfile = 'ldd_500m_v3SD.tif';
    accfile = 'accuflux_500m_v3SD_matlabGDAL.tif';
    mask500m=fullfile(root500m,"mask_union5km500m_v3SD"); % file extensions added later
    subbasinmaskfile='UIcatchment_union5km500m_cl_v3SD_matlabGDAL.tif';
    outletfile = 'outlets5_500m_v3SD.tif';
elseif scn=="orig"
    demfile = 'dem_500m.tif';
    dirfile = 'ldd_500m.tif';
    accfile = 'accuflux_500m.tif';
    mask500m=fullfile(root500m,"UIcatchment_union5km500m_mask");
    subbasinmaskfile='UIcatchment_union5km500m_cl_matlabGDAL.tif';
    outletfile = 'outlets5_500m.tif';
end

%% Filenames - 5km
root5km ="G:/SPHY/Sonu/full_runs_sanita/";

%% Create annual - QAll and Totr
if prepAnnualQ_T
    %% For annual - QAll
    Q_5km=fullfile(root5km, 'qall', 'QAllAvgY.map');
    convertmap2tif(Q_5km) %, Q_5kmY);
    cd(fullfile(root5km,'qall'))
    prefix="QAllAvgM";
    system(sprintf("pcrcalc %s13.map = (%s01.map + %s02.map + %s03.map + %s04.map + %s05.map + %s06.map + %s07.map + %s08.map + %s09.map + %s10.map + %s11.map + %s12.map)/12", prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix))
    system("pcrcalc Q13comp.map = QAllAvgY.map-QAllAvgM13.map")
    system("pcrcalc Q13pbias.map = (QAllAvgY.map-QAllAvgM13.map)/QAllAvgY.map*100")
    
    system("aguila Q13comp.map")
    convertmap2tif("QAllAvgM13.map") ;
    
    %% For annual - Totr
    cd(fullfile(root5km,'totr'))
    prefix="TotrAvgM";
    system(sprintf("pcrcalc %s13.map = (%s01.map + %s02.map + %s03.map + %s04.map + %s05.map + %s06.map + %s07.map + %s08.map + %s09.map + %s10.map + %s11.map + %s12.map)/12", prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix))
    convertmap2tif("TotrAvgM13.map") ;
    disp("Prepared annual Q and R files")
end

%% Load TotR and QAll .map, set nans as 0 and save as .mat
setnan=0;
if loadR_Q
    disp('Loading TotR and QAll');
    
    Rm_SPHY_mmday={};
    Qm_SPHY_m3sec={};
    for m = 1:13
        %% Convert TotR .map  files to .tiff
        R_5km=fullfile(root5km, 'totr',sprintf('TotrAvgM%02d.map',m));
        R_5kmtif=convertmap2tif(R_5km);
        Rm_SPHY_mmday{m} = loadSPHYtiff(R_5kmtif, nodata32bitFlo , showplot,setnan);
        
        %% Convert QAll .map  files to .tiff
        Q_5km=fullfile(root5km, 'qall', sprintf('QAllAvgM%02d.map',m));
        Q_5kmtif=convertmap2tif(Q_5km);
        Qm_SPHY_m3sec{m} = loadSPHYtiff(Q_5kmtif, nodata32bitFlo , showplot,setnan);
    end
    
    % Read and concatenate annual QAll fil too
    Q_5km=fullfile(root5km, 'qall', sprintf('QAllAvgY.map',m));
    Q_5kmtif=convertmap2tif(Q_5km);
    Qm_SPHY_m3sec{14} = loadSPHYtiff(Q_5kmtif, nodata32bitFlo , showplot,setnan);
    
    if exist('matfile')
        disp('Saving Rm and Qm data')
        matfile=fullfile(matpath,QRout);
        save(matfile,'-v7.3','Rm_SPHY_mmday', 'Qm_SPHY_m3sec');
        
    end
end


%% Load MASK file @ 500m w nans as 0
disp('Loading MASK and OUTLETs');

setnan=0;
basinmask = logical(loadSPHYtiff(strcat(mask500m, '_matlabGDAL.tif'), nodata8bit, showplot,setnan));
subbasinmask= loadSPHYtiff(fullfile(root500m, subbasinmaskfile), nodata32bitInt, showplot,setnan);
outside = ~basinmask;

% Basin outlets for catchment delineation, not outlet as nan
outlet = loadSPHYtiff(fullfile(root500m, outletfile), 0, showplot,setnan); % biggest subbasin in outlet3

%% Load acc @ 500m w nans as 1
disp('Loading ACC');
setnan=1;
acc = loadSPHYtiff(fullfile(root500m, accfile), nodata32bitFlo, showplot,setnan);

%% load SPHY DIR file @ 500m w nans as 5, i.e. SINK in SPHY
disp('Loading LDD');
setnan=5; % sink in sphy
lddSPHY = loadSPHYtiff(fullfile(root500m, dirfile), nodata8bit , showplot,setnan); % dir is uint8, cannot be nan

% Translate SPHY dir to Patrick dirs
fdir = translateSPHYdir2Patrick(lddSPHY);
adir = reversePatrickDirs(fdir);

%% Load dem @ 500m w nans as 0, no 0 in data
disp('Loading DEM');
setnan=0;
dem = loadSPHYtiff(fullfile(root500m, demfile), nodata32bitFlo, showplot,setnan);

%% Find channel
ma = nanmean(acc(:));
channel = acc>ma & basinmask==1; % only cells above mean flow accumulation area
showplot=1;
[rch,cch]=ind2sub(size(acc), find(channel));

%% Find drainage area for all 11 outlets
defdirs
showplot=0;
outID=unique(outlet); %~isnan(outlet)
catchment=0*acc;
c111=0;
r111=0;
ic111=0;
ir111=0;
for pt=outID(2:end)' %skip first as 0
    inlet=find(outlet==pt);
    basin = fastfindupstream_limSD(acc,fdir,drow,dcol, inlet,showplot,15000000);
    catchment = basin*pt+catchment;
    [ir0,ic0]=ind2sub(size(acc), inlet);
    ic111=[ic111; ic0];
    ir111=[ir111; ir0];
end

% Plot
catchment(basinmask==0)=nan;
figure
imagescnan(catchment);axis image; grid on
hold all; plot(ic111(2:end),ir111(2:end),"ok",'markersize',15)   %plot outlets
hold all; plot(cch(2:end),rch(2:end),'.b','markersize',1)   %plot channel
title('Subcatchments and channels in 500m')
legend('Outlets', 'All channels')
% Set cells outside basin to black color
c=colormap(parula);
c(1,:)=[0 0 0];
colormap(c)
%
disp("Catchments and mainstreams delineated")

%% Save files
if exist('matfile')
    disp('Saving acc, fdir, adir, dem and mask data')
    matfile=fullfile(matpath,data500m);
    save(matfile,'-v7.3','acc', 'fdir', 'adir','dem', 'basinmask', 'subbasinmask','outlet', 'outside','catchment');
end
disp("***************************************EOF***************************************")

