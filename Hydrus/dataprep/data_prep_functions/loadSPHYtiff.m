function data = loadSPHYtiff(fpath, nodataval , showplot,fillval)
% Read data in tiff and replace nodataval to a fillval (default = nan)
% Converts data to single type because of how imread handles tiff

if ~exist('fillval')
    fillval=nan;
end

% Get filename from path
[~, fname,fext] = fileparts(fpath);
data = single(imread(fpath)); % by default 16-bit floating-point images are returned as class single.

% fill nans and nodatavals w fillval
data(data==nodataval)=fillval;
data(isnan(data))=fillval;

if ndims(data)>2
    disp("Input tiff has multiple bands, cannot plot")
    return
end
%% Show data
if showplot
    figure;
    subplot(1,2,1)
    imagesc(data);axis image;
    colorbar
    
    grid on
    subplot(1,2,2)
    plot(sort(data(:)))
    grid on
    ylabel("Cell value")
    
    sgtitle(sprintf("Values in %s%s", fname,fext))
    fprintf("# of data cells in %s: %d \n", fname, sum(~isnan(data(:))))
end
end
