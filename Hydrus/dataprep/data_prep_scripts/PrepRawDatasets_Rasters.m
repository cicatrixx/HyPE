%% Load raw grid files. Project and clip them to the UI extents.
% Save into .mat formats.
% nodata values set to -9999 in tif files and 0 in .mat files
% For all global GCS grids,  reproject to pantpe, clip to UI extents, fill nodata=-9999 and
% save to tiff using gdal_wrap. Tiff is read using loadSPHYtiff() and saved
% to a mat file
clear 
close all

addpath(genpath(fullfile(pwd,'\Hydrus\devFiles')))
% run code from GitConnect fldr
root_data= fullfile(pwd,'\data\data_prep');
suffix="_UIprj"; % suffix added to output filename
matpath=fullfile(pwd,'\data\UI\data');
nodataval=-9999;
showplot=1;
setnan=0;
theorypotdata=0;

%% data relevant to theorypot
if theorypotdata
    %% 15s (500m) HydroSHEDs DEM in ESRI BIL format
    omat='Z500m_fill.mat';
    ifldr='as_dem_15s_bil';
    ifile=fullfile(root_data,ifldr,"as_dem_15s.bil");  %original downloaded from hydrosheds
    remap="near";
    gcs_wgs84= "EPSG:4326"; % the standard world GCS
    oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84); %% renamed to as_dem_500m_UIprj.tif and moved to the 500m folder
    
    % Resampled and reprojected DEM is filled in arcGIS and then loaded again
    %get Rw
    [~,Rw500] =geotiffread(oras);
    nodataval=32767; %16 bit signed
    omat='Z500m_fill.mat';
    ifldr='500m';
    ofile=fullfile(root_data,ifldr,"as_dem_500m_UIprj_fill.tif");
    data = loadSPHYtiff(ofile, nodataval , showplot, setnan);
    matfile=fullfile(matpath,omat);
    save(matfile,'data','Rw500')
    
    %% Merged 3s (100m) HydroSHEDs DEM grids for Indus
    omat='Z100m_fill.mat';
    ifldr='100m';
    ifile=fullfile(root_data,ifldr,"Indus_3s_con.tif");   % downloaded 3s bilcon from hydrosheds and merged them
    remap="near";
    gcs_wgs84= "EPSG:4326"; % the standard world GCS
    res=100;
    oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84,res);
    
    % Resampled and reprojected DEM is filled in arcGIS and then loaded
    ofile=fullfile(root_data,ifldr,"Indus_3s_con_UIprj_fill.tif");%get Rw
    [~,Rw100] =geotiffread(ofile);
    %
    nodataval=32767; %16 bit signed
    
    data = loadSPHYtiff(ofile, nodataval , showplot, setnan);
    %figure;imagesc(data);set(gca, 'ColorScale', 'log')
    matfile=fullfile(matpath,omat);
    save(matfile,'data','Rw100')
    
    disp("DEM 15s and 3s data loaded")
    
    %% ArcGIS based FDIR
    omat='FDIR500m.mat';
    ifldr='500m';
    ifile=fullfile(root_data,ifldr,"as_fdir_500m_UIprj02.tif");   % prepared in ArcGIS from hydrosheds DEM and altered in QGIS using Serval
    data = loadSPHYtiff(ifile, nodataval , showplot, setnan);
    matfile=fullfile(matpath,omat);
    save(matfile,'data')
    disp("FDIR data loaded")
    
    %% ArcGIS based ACC
    omat='ACC500m.mat';
    ifldr='500m';
    ifile=fullfile(root_data,ifldr,"as_facc_500m_UIprj02.tif");   % prepared in ArcGIS from hydrosheds DEM and altered in QGIS using Serval
    data = loadSPHYtiff(ifile, nodataval , showplot, setnan);
    matfile=fullfile(matpath,omat);
    save(matfile,'data')
    disp("acc data loaded")
    
    %% ArcGIS based outlets
    nodataval=-128; %8bitsigned
    omat='Outlets500m.mat';
    ifldr='500m';
    ifile=fullfile(root_data,ifldr,"outlets_500m_UIprj02.tif");   % prepared in ArcGIS from hydrosheds DEM and altered in QGIS using Serval
    data = loadSPHYtiff(ifile, nodataval , showplot, setnan);
    matfile=fullfile(matpath,omat);
    save(matfile,'data')
    disp("outlet data loaded")
    
    %% ArcGIS based catchments
    nodataval=255; %8bitunsigned
    omat='Catchments500m.mat';
    ifldr='500m';
    ifile=fullfile(root_data,ifldr,"catchments_500m_UIprj02.tif");   % prepared in ArcGIS from hydrosheds DEM and altered in QGIS using Serval
    data = loadSPHYtiff(ifile, nodataval , showplot, setnan);
    matfile=fullfile(matpath,omat);
    save(matfile,'data')
    disp("catchments data loaded")
end

%% ArcGIS based streamorder - STRAHLER
nodataval=-2147483648; %8bitunsigned
omat='streamorder_STRAHLER.mat';
ifldr='500m';
ifile=fullfile(root_data,ifldr,"streamorder_STRAHLER.tif");   % prepared in ArcGIS from hydrosheds DEM and altered in QGIS using Serval
data = loadSPHYtiff(ifile, nodataval , showplot, setnan);
matfile=fullfile(matpath,omat);
save(matfile,'data')
disp("streamorder data loaded")

%% ArcGIS based streamorder - STRAHLER Q based
nodataval=-2147483648; %8bitunsigned
omat='streamorder_STRAHLER_Qbased.mat';
ifldr='500m';
ifile=fullfile(root_data,ifldr,"streamorder_STRAHLER_Qbased.tif");   % prepared in ArcGIS from hydrosheds DEM and altered in QGIS using Serval
data = loadSPHYtiff(ifile, nodataval , showplot, setnan);
matfile=fullfile(matpath,omat);
save(matfile,'data')
disp("streamorder data loaded")

%% ArcGIS based streamorder - StreamOrder Q based
nodataval=127; %8bitsigned
omat='streamorder_Qclasses.mat';
ifldr='500m';
ifile=fullfile(root_data,ifldr,"streamorder_Qclasses_UIprj.tif");   % prepared in ArcGIS from hydrosheds DEM and altered in QGIS using Serval
data = loadSPHYtiff(ifile, nodataval , showplot, setnan);
matfile=fullfile(matpath,omat);
save(matfile,'data')
disp("streamorder data loaded")

%% ArcGIS based streamorder - SHREVE
nodataval=-2147483648; %8bitunsigned
omat='streamorder_SHREVE.mat';
ifldr='500m';
ifile=fullfile(root_data,ifldr,"streamorder_SHREVE.tif");   % prepared in ArcGIS from hydrosheds DEM and altered in QGIS using Serval
data = loadSPHYtiff(ifile, nodataval , showplot, setnan);
matfile=fullfile(matpath,omat);
save(matfile,'data')
disp("streamorder data loaded")

%% Save up georefrence
omat='pantpe_georef.mat';
pantpe500m=fullfile(root_data,"500m","extents_500m_v3SD.tif");
[~,Rw500m] =geotiffread(pantpe500m);
proj500m = geotiffinfo(pantpe500m);

pantpe100m=fullfile(root_data,"100m","Indus_3s_con_UIprj_fill.tif");
[~,Rw100m] =geotiffread(pantpe100m);
proj100m = geotiffinfo(pantpe100m);

matfile=fullfile(matpath,omat);
save(matfile,'Rw500m','proj500m','Rw100m','proj100m')

disp("Georef data loaded")

%% HYDE 3.2 for 2015 - .asc
%popd_<yr>.asc (population density, in inhabitants/km2 per gridcell)
nodataval=-9999; % added to function
omat='Popd.mat';
ifldr='HYDE_3_2_2015AD_pop';
ifile=fullfile(root_data,ifldr,"popd_2015AD.asc");
remap="near";
gcs_wgs84= "EPSG:4326"; % the standard world GCS
oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84);
%
data = loadSPHYtiff(oras, nodataval , showplot, setnan);
%figure;imagesc(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

disp("Pop data loaded")

%% Earthquake PGA grid - R2 use filled
nodataval=-3.402823466e+38;
omat='GSHAP.mat';
ifldr='gshap_globe_grd';
ifile=fullfile(root_data,ifldr,"gshap_globe_fill.tif");
remap="near";
gcs_wgs84= "EPSG:4326"; % the standard world GCS
oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84); %load pga
%
data = loadSPHYtiff(oras, nodataval , showplot);
%figure;imagesc(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

disp("Earthquake data loaded")

%% Landslide susceptibility grid
omat='LandSlideSusceptibility.mat';
ifldr='LHASA_LandSlideSusceptibilityMap';
ifile=fullfile(root_data,ifldr,"suscV1_1.tif");
remap="near";
gcs_wgs84= "EPSG:4326"; % the standard world GCS
oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84); %load pga
%
data = loadSPHYtiff(oras, nodataval , showplot, setnan);
%figure;imagesc(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

disp("Landslide susceptibility data loaded")

%% Landslide fuzzy probabilities - R2
omat='LandSlideSusceptibility_R2.mat';
ifldr='LHASA_UpdateFromStanley2021_CentralAsia';
ifile=fullfile(root_data,ifldr,"Unclassified.tif");
remap="near";
gcs_wgs84= "EPSG:4326"; % the standard world GCS
oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84); %load pga
%
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , showplot, setnan);
%figure;imagesc(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

disp("Landslide susceptibility data 2 loaded")

%% ESA Landuse grid - 2015
omat='Landuse.mat';
ifldr='ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7';
ifile=fullfile(root_data,ifldr,"ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif");
remap="near";
gcs_wgs84= "EPSG:4326"; % the standard world GCS
%
system(sprintf("gdalinfo %s",ifile))
oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84); %load pga
%
data = loadSPHYtiff(oras, nodataval , showplot, setnan);
%figure;imagesc(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

disp("Landuse data loaded")

% add LULC classes to fig
fid=fopen(fullfile(root_data,ifldr,"lccsclasses.txt"));
C=textscan(fid,'%d %s','Delimiter',';');
fclose(fid)
lclables=1;
figure(1);subplot(1,2,2)
set(gca,'YTick',C{1},'YTickLabel',C{2})

%% Tree density
omat='Tree density.mat';
ifldr='Crowther_Nature_Ecoregion';
ifile=fullfile(root_data,ifldr,"Crowther_Nature_Ecoregion.tif"); %
remap="near";
pcs= "+proj=igh"; % Goode Homolosine interrupted projected coordinate system
oras=grd2PantpeTiff(ifile,suffix,remap,pcs); %load in number of trees per km2 per cell
%
data = loadSPHYtiff(oras, nodataval , showplot, setnan);
%figure;imagesc(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

disp("Treecover data loaded")

%% Wouter-BIJL Water consumption data in mm/day reproject to 5km grids
omat='Dom_IndWC_mmday.mat';
ifldr='Wouter_WaterDemand_Bijl';
ifile=fullfile(root_data,ifldr,"totwc_1981-2010_mmday.nc"); %
remap="near";
gcs_wgs84= "EPSG:4326"; % the standard world GCS
res5km=5000; %5km in m
oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84,res5km);
%
data = loadSPHYtiff(oras, nodataval , 0,0);
%data = data*res5km*res5km; %convert to m3/day
figure;imagesc(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')
disp("Dom_Ind Water Consumption data loaded")

% %% Wouter-BIJL Water consumption data in m3/day convert nc to tiff...no reprojection
% ifldr='WaterDemand_Bijl_Wouter';
% ifile=fullfile(root_data,ifldr,"totwc_1981-2010_m3day.nc"); %
% oras=fullfile(root_data,ifldr,"totwc_1981-2010_m3day.tiff"); %
% system(sprintf("gdal_translate %s %s",ifile,oras))

%% Wouter- LPjML Irrigation Water demand in mm/day for 12 months reproject to 5km grids
omat='Irri_mmday.mat';
ifldr='Wouter_LPJML';
ifile=fullfile(root_data,ifldr,"irriwc_1981-2010_mmday.nc"); %
remap="near";
gcs_wgs84= "EPSG:4326"; % the standard world GCS
res5km=5000; %5km in m
oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84,res5km);

%
data = loadSPHYtiff(oras, nodataval , 0,0);
%data = data*res5km*res5km; %convert to m3/month
figure;imagesc(data(:,:,10));set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')
disp("LPjML Irrigation Water demand data loaded")

% Wouter- LPjML Irrigation Water demand in m3/day convert nc to tiff...no reprojection
ifldr='LPJML_Wouter';
ifile=fullfile(root_data,ifldr,"irriwc_1981-2010_m3day.nc"); %
oras=fullfile(root_data,ifldr,"irriwc_1981-2010_m3day.tiff"); %
system(sprintf("gdal_translate %s %s",ifile,oras))

%% Wouter- LPjML Irrigation Water CONSUMPTION in mm/day for 12 months reproject to 5km grids
omat='Irri_wc_mmday.mat';
ifldr='Wouter_Irrig_Consumption';
ifile=fullfile(root_data,ifldr,"irriwc_1981-2010_mmday_LTmonavg.nc"); %
remap="near";
gcs_wgs84= "EPSG:4326"; % the standard world GCS
res5km=5000; %5km in m
oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84,res5km);

%
data = loadSPHYtiff(oras, nodataval , 0,0);
%data = data*res5km*res5km; %convert to m3/month
figure;imagesc(data(:,:,10));set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')
disp("LPjML Irrigation Water CONSUMPTION data loaded")

% Wouter- LPjML Irrigation Water Consumption in m3/day convert nc to tiff...no reprojection
ifile=fullfile(root_data,ifldr,"irriwc_1981-2010_m3day_LTmonavg.nc"); %
oras=fullfile(root_data,ifldr,"irriwc_1981-2010_m3day_LTmonavg.tiff"); %
system(sprintf("gdal_translate %s %s",ifile,oras))


%% Wouter- LPjML Crop Yield data in g/m2/yr or tonne/km2-- absolute cell yield (crop_yield*crop_frac) BUT this output mat file is too large to open sometimes!
% orig data in gC/m2/yr, converted to g/m2/yr in R
omat='CropsPot_cell_g_m2_5minLPjML.mat';
ifldr='Wouter_LPJML_Yields';
remap="near";
gcs_wgs84= "EPSG:4326"; % the standard world GCS
showplot=0;

% Loop through files for 32 crop types
for ctype=4:32
    ifile=fullfile(root_data,ifldr,sprintf("yield_g_m2_yr_TS_1981_2010_pft_%d.nc",ctype)); %
    oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84);
    % load tiff and convert nodataval to 0
    Crops_data{ctype} =  loadSPHYtiff(oras, nodataval , showplot);     % potential production tons(wet)/km^2
    %figure;imagescnan(Crops_data{ctype}(:,:,1)); set(gca, 'ColorScale', 'log')
end
matfile=fullfile(matpath,omat);
save(matfile,'Crops_data', '-v7.3')

disp("New crop potentials data loaded")

%% IGNORE Wouter- LPjML Crop Yield data in g/m2/yr or tonne/km2 -- harvest (crop_yield)
% orig data in gC/m2/yr, converted to g/m2/yr in R
omat='CropsPot_harvest_g_m2_5minLPjML.mat';
ifldr='LPJML_Wouter_Yields';
remap="near";
gcs_wgs84= "EPSG:4326"; % the standard world GCS
showplot=0;

% Loop through files for 32 crop types
for ctype=2:32
    ifile=fullfile(root_data,ifldr,sprintf("harvest_g_m2_yr_TS_1981_2010_pft_%d.nc",ctype)); %
    oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84);
    % load tiff and convert nodataval to 0
    Crops_data{ctype} =  loadSPHYtiff(oras, nodataval , showplot);     % potential production tons(wet)/km^2
    %figure;imagescnan(Crops_data{ctype}(:,:,1)); set(gca, 'ColorScale', 'log')
end
matfile=fullfile(matpath,omat);
save(matfile,'Crops_data', '-v7.3')

disp("New crop potentials--harvest based data loaded")

%% OLD Crop potentials .asc
omat='CropsPot.mat';
ifldr='PotentialYield_IMAGE';
remap="near";
gcs_wgs84= "EPSG:4326"; % the standard world GCS

for i=1:18 % Only rainfed food crops
    ifile=fullfile(root_data,ifldr,sprintf('CropPot_actual_potential_yield_2011_%d.asc',i));
    oras=grd2PantpeTiff(ifile,suffix,remap,gcs_wgs84);
    % load tiff and convert nodataval to 0
    Crops_data{i} =  loadSPHYtiff(oras, nodataval , showplot, setnan);     % potential production tons(wet)/km^2
    
end
%figure;imagesc(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'Crops_data')

disp("Crop potentials data loaded")
