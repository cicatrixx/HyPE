%% Compile theoretical pot related files into matlab arrays
% based on processSPHYoutputs but this one is for ArcGIS based outputs and Qall
% Compiled rasters already preloaded into matlab
% fdir is changed from ESRI to Patrick dir metrics
% Things are kept in original extents for now. No masking of endorheic and
% small basins

% Created By    : Sanita Dhaubanjar on 08 Feb 2021
% Created For	: SustaIndus project WP2
%=========================
clc
close all
clear all

root500m=fullfile(pwd,'data','UI','data');
dataout500m='UI500m_ArcGIS.mat'; % for Z, acc, fdir and adir
matout=fullfile(root500m,dataout500m);

%% Filenames - 500m
demfile = 'Z500m_fill.mat';
dirfile = 'FDIR500m.mat';
accfile = 'ACC500m.mat';
outletfile = 'Outlets500m.mat';
catchmentfile='Catchments500m.mat';

load(fullfile(root500m,demfile),'data')
dem=data;
load(fullfile(root500m,dirfile),'data')
dir=data;
load(fullfile(root500m,accfile),'data')
acc=data+1; %the new acc from ArcGIS has cellval=0. So modified acc to be acc+1 to avoid any potential errors in david?s other codes that assume acc>0
load(fullfile(root500m,outletfile),'data')
outlets=data;
load(fullfile(root500m,catchmentfile),'data')
catchments_GIS=data;

%% Create mask as out catchment
basinmask=catchments_GIS~=0;
outside=~basinmask;

%% Create fdir and adir in patrick direction metrics
% Translate ESRI dir to Patrick dirs
dir(outside)=247;
fdir = translateESRIdir2Patrick(dir);
adir = reversePatrickDirs(fdir);

%% Check if DEM has any negative heads
[~, cell_Hgross]=evalCell2CellTheoreticalPot(dem, fdir, dem*0.1, 0);
if any(cell_Hgross(:)<0)
    disp('Check dem for negative heads!!')
    keyboard
end

%% Create catchment area for all 11 outlets
run('defdirs.m')
showplot=0;
outID=unique(outlets); %~isnan(outlet)
catchments=0*acc;
c111=0;
r111=0;
ic111=0;
ir111=0;
for pt=outID(2:end)' %skip first as 0
    inlet=find(outlets==pt);
    basin = fastfindupstream_limSD(acc,fdir,drow,dcol, inlet,showplot,15000000);
    catchments = basin*pt+catchments;
    [ir0,ic0]=ind2sub(size(acc), inlet);
    ic111=[ic111; ic0];
    ir111=[ir111; ir0];
end

%% Create channel - assuming all basin boundary!
[rch,cch,~,channel]=getChannel(acc, catchments==0);
disp('Created catchments and channels')

%% Plot
catchments(basinmask==0)=0;
figure
imagescnan(catchments);
hold all; plot(ic111(2:end),ir111(2:end),"ok",'markersize',15)   %plot outlets
hold all; plot(cch(2:end),rch(2:end),'.w','markersize',1)   %plot channel
%title('New Catchment')
title('New subcatchments and channels in 500m')
legend('Outlets', 'All channels')

disp("Catchments and mainstreams delineated")

if all(all(catchments==catchments_GIS))
    disp("Matlab and arcGIS catchments are the same")
else
    disp("Check catchment delineation!")
    keyboard
end

%% Rename basins and outlets
cleanupIDs=1;
if cleanupIDs
    % Rename outlet
    oldIDs=[0,1:5,8,10,12,13,14,15,16];
    newIDs=[0,102,101,103,109,110,104,105,111,106, 107,108,112];
    %
    catchments2=changem(catchments,newIDs,oldIDs);
    outlets2=changem(outlets,newIDs,oldIDs);
    
    %for saving - unique IDs and names
    basinIDs=[101:112]';
    basinnames={'Kabul','Swat','Indus Main','Jhelum','Chenab','Ravi','Beas','Satluj','Misc','Misc','Misc','Endorheic'}';
    outlet_idx=find(outlets2>0);% & outlets2 ~>108); %find only outlets that matter
    ptID=outlets2(outlet_idx);
    %order the indexes based on basinID order
    outlet_idxs=changem(basinIDs,outlet_idx,ptID);
    %order old ids based on basinID order
    old_basinIDs=changem(basinIDs,oldIDs,newIDs);
    basinlabels=table(basinIDs,basinnames,outlet_idxs,old_basinIDs);
    
    %% Plot catchments and outlets w IDs as well as name for double checking
    % create catchment map and colors for it
    cmap = linspecer(length(basinIDs),'qualitative'); %
    
    %get basin names in the same order as ptID
    %basin4ptID=categorical(ptID,basinIDs,basinnames);
    
    [ir0,ic0]=ind2sub(size(outlets2), outlet_idxs);
    
    f=figure('Position', get(0, 'Screensize'));
    imagescnan(maskBasin(catchments2, basinmask))
    text(ic0, ir0,num2str(basinIDs),'FontSize' ,8,'BackgroundColor',.9*[1 1 1]);%,'FontWeight','bold'
    text(ic0, ir0,basinnames,'FontSize' ,8,'BackgroundColor',.9*[1 1 1]);%,'FontWeight','bold'
    hold all; plot(ic0,ir0,".k",'markersize',15)   %plot outlets
    colormap(cmap)
    colorbar('Ticks', basinIDs,'TickLabels',basinnames,'Direction','reverse')
    grid on
    title("Outlet IDs and Catchments")
    %orient('landscape')
    %print('G:\GitConnect\output\TheoryPot_WrapUp\UnionCatchments_Named.pdf','-dpdf','-fillpage')
    
    %rename
    outlet_orig=outlets;
    outlets=outlets2;
    catchments=catchments2;
end

%% Create flowdist
flowdist=flowdistanceSD(acc,fdir,outside);
figure; imagescnan(flowdist)
title('New flowdist')
%% Compare w old
if showplot
    old=load('G:\SurfDrive\GitConnect\data\UI\data\UI500m.mat', 'fdir','acc','basinmask_union','catchment');
    old.fdist=flowdistanceSD(old.acc,old.fdir, ~old.basinmask_union);
    figure; imagescnan(old.fdist)
    title('Old flowdist')
    
    figure; imagescnan(old.catchment)
    title('Old catchment')
    ma = nanmean(old.acc(:));
    old.channel = old.acc>ma & old.basinmask_union==1; % only cells above mean flow accumulation area
    [old.rch,old.cch]=ind2sub(size(old.acc), find(old.channel));
    hold all; plot(old.cch(2:end),old.rch(2:end),'.b','markersize',1)   %plot channel
end
%% Save files
if exist('matout')
    disp('Saving acc, fdir, adir, dem, basinnames and mask data')
    save(matout,'-v7.3','acc', 'fdir', 'adir','dem','flowdist', ...
        'basinmask','catchments', 'basinlabels', 'outlets', 'outlet_orig');% all based on renamed versions
    disp('Processed data is saved')
end

disp("***************************************EOF***************************************")