%% Load raw shapefiles. Project, clip and rasterize them to the UI extents. 
% Save into .mat formats.
% nodata values set to -9999 in tif files and 0 in .mat files
clear all
close all
%root_data='D:\GitConnect\data\data_prep';
root_data='G:\SurfDrive\GitConnect\data\data_prep';
suffix="_UIprj"; % suffix added to output filename
matpath='G:/SurfDrive/GitConnect/data/UI/data';

%% Country boundaries
omat='Countries.mat';
ifldr='GADM_Countries';
ifile=fullfile(root_data,ifldr,"Indus_Countries.shp");
selattribute ='ISO_ID' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
oras=shp2tiff(ifile,selattribute,suffix);
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , 1,0);

%figure;imagesc(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% Protected Areas
omat='WDPA.mat';
ifldr='WDPA_Jun2020';
ifile=fullfile(root_data,ifldr,"WDPA_Jun2020-shapefile-polygons.shp");
selattribute ='WDPAID' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
oras=shp2tiff(ifile,selattribute,suffix);
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , 1,0);

%figure;imagesc(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')
    
%% Water Bodies - Natural Only
omat='WaterBodies.mat';
ifldr='HydroLAKES_polys_v10_shp';
ifile=fullfile(root_data,ifldr,"HydroLAKES_polys_v10_NaturalOnly.shp");
selattribute ='Hylak_id' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
oras=shp2tiff(ifile,selattribute,suffix);
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , 1,0);

%figure;imagesc(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% Faults and Thrusts
omat='Seismic_Fault_Thrust.mat';
ifldr='Central_Asia_Fault_Database';
ifile=fullfile(root_data,ifldr,"cafd_faults_27072015_GCS.shp");
selattribute ='ID' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
oras=shp2tiff(ifile,selattribute,suffix);
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , 1,0);

%figure;imagesc(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% World Powerlines - replaced!
% omat='Powerlines.mat';
% ifldr='World_powerlines';
% ifile=fullfile(root_data,ifldr,"World_powerlines.shp");
% selattribute ='id' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
% oras=shp2tiff(ifile,selattribute,suffix);
% nodataval=-9999;
% data = loadSPHYtiff(oras, nodataval , 1,0);
% 
% %figure;imagescnan(data);set(gca, 'ColorScale', 'log')
% matfile=fullfile(matpath,omat);
% save(matfile,'data')

%% TL_fromOverpassTurbo
omat='Powerlines.mat';
ifldr='TL_fromOverpassTurbo';
ifile=fullfile(root_data,ifldr,"TransmissionLine_UI2.shp");
selattribute ='OID_' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
oras=shp2tiff(ifile,selattribute,suffix);
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , 1,0);

%figure;imagescnan(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% Mountain areas
omat='MountainAreas.mat';
ifldr='GMBAv1.2';
ifile=fullfile(root_data,ifldr,"GMBAMountainInventory_v1_2-Asia.shp");
selattribute ='id' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
oras=shp2tiff(ifile,selattribute,suffix);
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , 1,0);

%figure;imagescnan(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% Glaciers
omat='Glaciers.mat';
ifldr='RGIv60_SouthAsiaWest';
ifile=fullfile(root_data,ifldr,"14_rgi60_SouthAsiaWest.shp");
selattribute ='id' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
oras=shp2tiff(ifile,selattribute,suffix);
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , 1,0);

%figure;imagescnan(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% Glacial lakes
omat='GlacialLakes.mat';
ifldr='ICIMOD_HKH_GlacialLake_2005';
ifile=fullfile(root_data,ifldr,"GlacialLake_5basins_HKH.shp");
selattribute ='ID' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
oras=shp2tiff(ifile,selattribute,suffix);
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , 1,0);

%figure;imagescnan(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% Roads - GRIP4
omat='Roads.mat';
ifldr='GRIP4';
ifile=fullfile(root_data,ifldr,"GRIP4_region6.shp");
selattribute ='select4LI' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
%oras=shp2tiff(ifile,selattribute,suffix);
%nodataval=-9999;

%% Roads - GRIP4 Use edited road map for convert to raster
omat='Roads2.mat';
ifile=fullfile(root_data,ifldr,"GRIP4_region6_UIprj_edited4.shp");
nodataval=-9999;

% convert to tiff
% For rasterization
xres=500;
yres=500;
nodataval=-9999;

% Pantpe bounds for UI
xmin=-1670000.0;
xmax=-230000.0;
ymin=30000.0;
ymax=840000.0;
pantpe_bbox=num2str([xmin, ymin, xmax, ymax]);
ofileras=strrep(ifile,".shp",strcat(suffix,".tiff"));
system(strjoin(["gdal_rasterize -a",selattribute,"-tr",xres,yres,"-te",pantpe_bbox, "-a_nodata", nodataval, ifile, ofileras]));

% save to mat
data = loadSPHYtiff(ofileras, nodataval , 1,0);
figure;imagescnan(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% GRanD
omat='Grand.mat';
ifldr='Grand';
ifile=fullfile(root_data,ifldr,"grand_dams_David.shp");
selattribute ='GRAND_ID' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
oras=shp2tiff(ifile,selattribute,suffix);
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , 1,0);

%figure;imagescnan(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% Settlements
omat='Settlements.mat';
ifldr='Settlements_ICIMOD';
ifile=fullfile(root_data,ifldr,"populated_place.shp");
selattribute ='ID' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
oras=shp2tiff(ifile,selattribute,suffix);
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , 1,0);

figure;imagescnan(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% Cultural heritages
omat='CulturalHeritage.mat';
ifldr='HeritageSites';
ifile=fullfile(root_data,ifldr,"HeritageSites_compiled.shp");
selattribute ='OBJECTID' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
oras=shp2tiff(ifile,selattribute,suffix);
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , 1,0);

figure;imagescnan(data);set(gca, 'ColorScale', 'log')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% GLOF Hazards
omat='potentialGLOFs.mat';
ifldr='GLOFrisk';
ifile=fullfile(root_data,ifldr,"GLOF_gl_Indus_2001_2.shp");
selattribute ='SN' ; % rasterize based on this val column, here Unique identifier for parcels or zones within a protected area.
oras=shp2tiff(ifile,selattribute,suffix);
nodataval=-9999;
data = loadSPHYtiff(oras, nodataval , 1,0);

figure;imagescnan(data);
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% GLOF Hazards Use edited points for convert to raster
omat='potentialGLOFs2.mat';
ifldr='GLOFrisk';
selattribute ='SN' ; % rasterize based on this val column
ifile=fullfile(root_data,ifldr,"GLOF_gl_Indus_2001_2_UIprj2.shp");
nodataval=-9999;

% convert to tiff
% For rasterization
xres=500;
yres=500;
nodataval=-9999;

% Pantpe bounds for UI
xmin=-1670000.0;
xmax=-230000.0;
ymin=30000.0;
ymax=840000.0;
pantpe_bbox=num2str([xmin, ymin, xmax, ymax]);
ofileras=strrep(ifile,".shp",strcat(suffix,".tiff"));
system(strjoin(["gdal_rasterize -a",selattribute,"-tr",xres,yres,"-te",pantpe_bbox, "-a_nodata", nodataval, ifile, ofileras]));

% save to mat
data = loadSPHYtiff(ofileras, nodataval , 0,0);
figure;imagescnan(data);set(gca, 'ColorScale', 'linear')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%% HP database Use projected points for convert to raster
omat='visualizedHP.mat';
ifldr='500m';
selattribute ='PourID' ; % rasterize based on this val column
ifile=fullfile(root_data,ifldr,"HP_forMLAB_UIprj.shp");
nodataval=-9999;

% convert to tiff
% For rasterization
xres=500;
yres=500;
nodataval=-9999;

% Pantpe bounds for UI
xmin=-1670000.0;
xmax=-230000.0;
ymin=30000.0;
ymax=840000.0;
pantpe_bbox=num2str([xmin, ymin, xmax, ymax]);
ofileras=strrep(ifile,".shp",strcat(suffix,".tiff"));
system(strjoin(["gdal_rasterize -a",selattribute,"-tr",xres,yres,"-te",pantpe_bbox, "-a_nodata", nodataval, ifile, ofileras]));

% save to mat
data = loadSPHYtiff(ofileras, nodataval , 0,0);
figure;imagescnan(data);set(gca, 'ColorScale', 'linear')
matfile=fullfile(matpath,omat);
save(matfile,'data')

%%
function ofileras=shp2tiff(ifile,selattribute,suffix,iproj4)
% Project input shapefile to pantpe projection and clip to all Indus
% bounding box. Rasterize projected shapefile and clip to Upper Indus box. 
% Output raster is 500m .tiff w suffix added to ifilename.
% Uses GDAL tools

% Proj4 string for pantpe
pantpe_proj4='"+proj=aea +lat_1=26 +lat_2=38 +lat_0=30 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"';

% For full Indus w some lee way
lonmin=56;
lonmax=90;
latmin=20;
latmax=42;
indus_bbox=num2str([lonmin, latmin, lonmax, latmax]);

% Pantpe bounds for UI
xmin=-1670000.0;
xmax=-230000.0;
ymin=30000.0;
ymax=840000.0;
pantpe_bbox=num2str([xmin, ymin, xmax, ymax]);

% For rasterization
xres=500;
yres=500;
nodataval=-9999;

ofileprj=strrep(ifile,".shp",strcat(suffix,".shp"));
ofileras=strrep(ifile,".shp",strcat(suffix,".tiff"));
%% Project shp file to pantpe
%system(strjoin('gdalinfo ',ifile))
%using spat instead of clipsrc selected all features in .shp that intersect
%with or fall within the bounding box
system(strjoin(["ogr2ogr -spat",indus_bbox,"-t_srs",pantpe_proj4, ofileprj, ifile]))

%[status,result] = ;assert(status>0,result)
disp("Projected and clipped shapefile")

%% Rasterize shp file, selecting a specific attribute to convert from
% Update pixels whose center point is within the polygons in shp file
system(strjoin(["gdal_rasterize -a",selattribute,"-tr",xres,yres,"-te",pantpe_bbox, "-a_nodata", nodataval, ofileprj, ofileras]));
disp("Rasterized and clipped shapefile")
end