function savemat2Pantpetiff(tiffoutpath, matrix_data, fpath2Rw)
% saves matrix in Pantpe PCS projection specified by Rw and proj as geotiff
% Manually set nan vals to -9999 as not sure how matlab handles this otherwise.

% Load Indus Rw path if not specified by user
if ~exist('fpath2Rw','var')
    fpath2Rw=fullfile(pwd,'\data\ASIA\Basin_UIB\PantpeBasin_3.mat');
end
    setnodatavals=-9999;
    matrix_data(isnan(matrix_data))=setnodatavals;
    load(fpath2Rw,'Rw','proj')
    geotiffwrite(tiffoutpath,matrix_data,Rw,'GeoKeyDirectoryTag',proj.GeoTIFFTags.GeoKeyDirectoryTag)

    %Sometimes Arc doesn't recognize this as the NoData value. To fix this, 
    % before adding the GeoTiff as a layer in your map, 
    % right-click the layer in the Catalog, click properties, and then 
    % look for 'NoData Value' under 'Raster Information'. 
    % Click the 'Edit...' button, and then change the NoData Value to -9999.
end