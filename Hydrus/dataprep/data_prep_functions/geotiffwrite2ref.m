function geotiffwrite2ref(ofname, idata, georefname)
% write idata to a geotiff file assuming the georeferencing in reference file.
% georefname and ofname are fullpaths to reference file and output geotiff.

% Get georef information
[~,R] =geotiffread(georefname);
Rinfo = geotiffinfo(georefname); 

% Write data to geotiff
geotiffwrite(ofname,idata,R,'GeoKeyDirectoryTag',Rinfo.GeoTIFFTags.GeoKeyDirectoryTag)

end