function maskdata=maskBasin(idata, mask, fillval)
% Return matrix same size as data with filled value in all cells that have mask=0
% Nans in idata as not changed to fillval! so better to remove nans in
% idata
if ~exist('mask','var')
    fname='G:\PaperData\DavidGernaat\Sanita_model_package\data\ASIA\Basin\IndusAll.mat';
    load(fname,'indusall');
    mask=indusall;
end

%% Replace cells where mask=0
maskdata=idata;
for i=1:size(idata,3)
    tmp=idata(:,:,i);

    if exist('fillval','var')
        tmp(mask==0)=fillval;
    else
        tmp(mask==0)=nan;
    end
    maskdata(:,:,i)=tmp;
end
