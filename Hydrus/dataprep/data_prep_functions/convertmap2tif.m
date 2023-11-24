function ofname=convertmap2tif(ifname, ofname)
%% Convert a .map R files to .tiff using gdal_translate
% ifname and ofname are fullpaths to the input and output files
% If ofname is unspecified, ifname is suffixed w _matlabGDAL to get ofname

if ~exist("ofname")
    ofname = strrep(ifname,".map","_matlabGDAL.tif");  
end

system(sprintf("gdal_translate -of GTiff %s %s", ifname,ofname));
end
