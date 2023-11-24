function h = imagescnan( img_data, bgcolor )
% a wrapper for imagesc, with some formatting going on for nans
% plotting data. Removing and scaling axes (this is for image plotting)
% FROM:
% https://stackoverflow.com/questions/14933039/matlab-imagesc-plotting-nan-values
% and https://nl.mathworks.com/matlabcentral/answers/73492-assigning-different-color-to-nan-values-in-2d-matrix

if ~exist('bgcolor','var')
     % third parameter does not exist, so default it to something-- but using nargin is much much faster than using exist.
      bgcolor = 1*[1 1 1];
 end

h = imagesc(img_data);
axis image
%turn off X and Y axis ticks and labels
set(gca,'XTick',[],'YTick',[])

% setting axis background color
set(gca, 'color', bgcolor)

% setting alpha values as 0(=transparent) for nans and 1(=opaque) for nonnans
if ismatrix( img_data )
  set(h, 'AlphaData', ~isnan(img_data))
elseif ndims( img_data ) == 3
  set(h, 'AlphaData', ~isnan(img_data(:, :, 1)))
end

if nargout < 1
  clear h
end

