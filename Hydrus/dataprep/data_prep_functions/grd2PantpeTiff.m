function ofile=grd2PantpeTiff(ifile,suffix,remap,src_prj,res,t_prj)
% Project input raster to pantpe projection and clip to Upper Indus box.
% Target projection is pantpe if not specified otherwise. Nodataval for out file set as -9999.
% Output raster is 500m or user-specified res, .tiff w suffix added to ifilename.
% Uses GDALWARP: https://gdal.org/programs/gdalwarp.html

% For rasterization
nodataval=-9999; % nodataval for output file
if ~exist('res')
    res = 500;
end

% Pantpe bounds for UI
xmin=-1670000.0;
xmax=-230000.0;
ymin=30000.0;
ymax=840000.0;
pantpe_bbox=num2str([xmin, ymin, xmax, ymax]);
% Proj4 string for pantpe
pantpe_proj4='"+proj=aea +lat_1=26 +lat_2=38 +lat_0=30 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"';

if ~exist('t_prj'); t_prj=pantpe_proj4; end

% Get filename from path
[fpath, fname,fext] = fileparts(ifile);
ofile=fullfile(fpath,strcat(fname,suffix,".tif"));

%% Project and clip grd to pantpe
%gdalwarp -s_srs EPSG:4326 -t_srs -tr 500.0 500.0 -r near -te -1670000.0 30000.0 -230000.0 840000.0 -of GTiff popd_2015AD.asc reprj_crop_mlab.tif
system(strjoin(["gdalwarp","-s_srs", src_prj, "-t_srs",t_prj,"-tr",res,res,"-r",remap,"-te",pantpe_bbox,"-dstnodata", nodataval, "-of GTiff", ifile, ofile]));
disp("Projected and clipped grid")
end