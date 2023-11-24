%% Modified David's downscaleQ_NAM_loop to use preloaded 5km and 500m R to get Q
% downscaleQ_wWC.m
% Created By    : Sanita Dhaubanjar on 24 Nov 2020
% Created For	: SustaIndus project WP2
% Last Updated  : 22 Feb 2022

% uses low res runoff, irrigation demand and domestic+industrial water consumption
% in mm/day to generate Q and Qwc in m3/day

clc
clear all
close all
runtests=0;
root=pwd;
matpath=fullfile(root,"data","UI","data");
saveQ2mat=fullfile(matpath,"MLAB500m_40yrClimatology_m3day_R4.mat");
createplots=1;
addpath(genpath("G:\SurfDrive\GitConnect\Hydrus"),"G:\SurfDrive\GitConnect\Hydrus\devFiles\")
%% Load 500m basin data acc, adir
disp('Loading acc and dir at 500m');
data500m='UI500m_ArcGIS.mat';
load(fullfile(matpath,data500m), 'acc', 'adir','outside','basinmask','channel') %,  'fdir', 'outlet'
% nan areas outside the basin to ignore them
%acc(outside)=nan;

%% Load preloaded 5km runoff maps
QRout='SPHY5km_40yrClimatology.mat';
load(fullfile(matpath,QRout),'Rm_SPHY_mmday');
%Rm=Rm_SPHY_mmday;
% Convert R at 5km from mm/day to m3/day
for m=1:13
    R5kmmat_m3day{m} = Rm_SPHY_mmday{m}*5000^2*1e-3; %% from mm/day to m3/day at 5000m or 5km
end

%% Load 5km Water consumption hil - same for all months!
disp('Loading wc hil');
load(fullfile(matpath,'Dom_IndWC_mmday.mat'));
mhilmmday=data;
% Load 5km Water consumption irr
disp('Loading wc irr');
load(fullfile(matpath,'Irri_wc_mmday.mat'));
irrmmday=data;

%Eval annual avg
irrmmday(:,:,13) = sum(data(:,:,1:12),3)/12;

% Sum 5km irri and hil and convert nans to 0
for m=1:13
    tmp = irrmmday(:,:,m) + mhilmmday;
    tmp(isnan(tmp))=0;
    WCt_mmday{m}=tmp;
end

%% get Q natural
scale_mm_to_m3day = 500^2*1e-3;  %for 500m cell size
[Q500m_m3day, R500m_m3day]=downscaleR2Q(Rm_SPHY_mmday, acc, adir, scale_mm_to_m3day);
disp('Routed Qnatural');

%% get Qwc
[Qwc_500m_m3day, Rwc_500m_m3day]=downscaleR2Q(Rm_SPHY_mmday, acc, adir, scale_mm_to_m3day,WCt_mmday);
disp('Routed Qwc');

%% Save mat file
save(saveQ2mat,'-v7.3','Q500m_m3day', 'R500m_m3day', 'Qwc_500m_m3day', 'Rwc_500m_m3day');
disp('Qnatural and Qwc files saved!');

 %% Saved annual Q and R to geotiff
    savemat2Pantpetiff(fullfile(root,"data","data_prep","500m","Q500m_M13_m3s_R4.tif"), Q500m_m3day{13}/(24*60*60))
    savemat2Pantpetiff(fullfile(root,"data","data_prep","500m","R500m_M13_m3s_R4.tif"), R500m_m3day{13}/(24*60*60))
    disp('Annual Qnatural and Qwc files saved as tiff!');

%% Compare previous Q generation w now after adding wc codes
if runtests
    m=13;
    oldQ=load('G:\SurfDrive\GitConnect\data\UI\data\MLAB500m_40yrClimatology_m3day.mat');
    %
    [newwWCcode_Q500m_m3day, newwWCcode_R500m_m3day]=downscaleR2Q(Rm_SPHY_mmday, acc, adir, 500^2*1e-3);
    all(oldQ.Q500m_m3day{m}==newwWCcode_Q500m_m3day{m},'all')
    all(oldQ.R500m_m3day{m}==newwWCcode_R500m_m3day{m},'all')
    disp('Compared old and new discharge')
    
    %% Rerun w dummy WC - 0 WC and WC=30% of R
    for m=1:13
        Wcct{m}=0*Rm{m};
    end
    [dWC_Q500m_m3day, dWC_R500m_m3day]=downscaleR2Q(Rm_SPHY_mmday, acc, adir, 500^2*1e-3,Wcct);
    all(oldQ.Q500m_m3day{m}==dWC_Q500m_m3day{m},'all')
    
    %
    for m=1:13
        Wcct{m}=0.3*Rm{m};
    end
    [dWC_Q500m_m3day, dWC_R500m_m3day]=downscaleR2Q(Rm_SPHY_mmday, acc, adir, 500^2*1e-3,Wcct);
    
    % all cells have 30% reduction in water consumption
    Qdiff=(oldQ.Q500m_m3day{m}-dWC_Q500m_m3day{m})./oldQ.Q500m_m3day{m}*100;
    Qdiff(~basinmask)=0;
    figure;imagesc(Qdiff)
end


%% Plot: compare Qnatural and Qwc for annual

if createplots
    m=13;
    Qdiff=(Q500m_m3day{m}-Qwc_500m_m3day{m})./Q500m_m3day{m}*100;
    
    figure;
    subplot(2,3,1)
    imagescnan(maskBasin(Q500m_m3day{m},basinmask))
    title('Qnatural in m^3/day')
    
    subplot(2,3,2)
    imagescnan(maskBasin(Qwc_500m_m3day{m},basinmask))
    title('Qwc in m^3/day')
    
    subplot(2,3,3)
    imagescnan(maskBasin(Qdiff,basinmask))
    title('Percentage reduction in Q after WC')%'(Qnatural-Qwc)/Qnatural*100')
    
    subplot(2,3,5)
    imagescnan(WCt_mmday{m}) % cant mask this as this is in 5km!
    set(gca, 'ColorScale', 'log')
    title('Total irri dom ind demand in mm/day')
    sgtitle("For m=13")
end

%% Plot: with all 13 maps - Q generally bad but good for quick overview
defdirs
if createplots
    nTS = length(Rm_SPHY_mmday);
    Qmax = max(cell2mat(Q500m_m3day),[],'all');
    figure
    mname = month(datetime(2013,01,01):calmonths(1):datetime(2013,12,31),'name');
    
    for i=1:12
        subplot(3, 6,6 + i)
        imagescnan(maskBasin(Qwc_500m_m3day{i},basinmask))
        grid on
        caxis([0, Qmax])
        title(mname{i})
        set(gca, 'ColorScale', 'log')
    end
    subplot(3, 6,3:4)
    imagescnan(maskBasin(Qwc_500m_m3day{13},basinmask))
    grid on
    caxis([0, Qmax])
    c = colorbar('eastoutside');
    c.Label.String = 'Discharge (m^3/day)';
    set(gca, 'ColorScale', 'log')
    title("Annual Average")
    sgtitle('Q in m^3/day')
end

%% Find drainage area and mainstream for all 11 outlets
if createplots
    Qin=Q500m_m3day{13}; 
    %% Find channel
%     ma = nanmean(acc(:));
%     channel = acc>ma;% & basinmask==1; % only cells above mean flow accumulation area
    [rch,cch]=ind2sub(size(acc), find(channel));
  
    %     clipQ=nan*Qwc;
    %     clipQ(500:1000,2000:2500)=Qwc(500:1000,2000:2500);
    %     chkbufbasin=createBuffer(selbasin>0,50);
    %     tmpQ=Qwc;
    %     tmpQ(~chkbufbasin)=nan;
    load(fullfile(matpath,data500m),  'fdir', 'outlets', 'basinmask')
    showplot=0;
    outID= unique(outlets(outlets(:)~=0));%&outlets(:)~=110));
    catchments1=0*acc;
    c111=0;
    r111=0;
    ic111=0;
    ir111=0;
    
    % ID indices for mainstream cells based on current discharge
    for pt=outID'
        inlet=find(outlets==pt);
        basin = fastfindupstream_limSD(acc,fdir,drow,dcol, inlet,showplot,15000000);
        catchments1 = basin*pt+catchments1;
        
        mainstream = find_mainstream(Qin,inlet,adir);
        [r0,c0]=ind2sub(size(acc), mainstream);
        [ir0,ic0]=ind2sub(size(acc), inlet);
        c111=[c111; c0];
        r111=[r111; r0];
        ic111=[ic111; ic0];
        ir111=[ir111; ir0];
    end
    
    
    % Plot
    % catchments1(basinmask==0)=nan;
    %catchments1(catchments1==-999)=-10;
    figure
    imagescnan(maskBasin(catchments1,basinmask));axis image; grid on
    hold all; plot(c111(2:end),r111(2:end),'.r','markersize',12)  %plot mainstream
    hold all; plot(ic111(2:end),ir111(2:end),"ok",'markersize',15)   %plot outlets
    hold all; plot(cch(2:end),rch(2:end),'.b','markersize',1)   %plot channel
    title('Subcatchments and channel in 500m')
    legend('Mainstream - highest Q based', 'Outlets', 'All channels')
    
    disp("Catchments and mainstreams delineated")
end
%% Find and plot mainstream - one outlet only
% if createplots
%         %[~,inlet]=max(acc(:)); % cell w highest acc
%     inlet=find(outlets==103); % biggest subbasin in outlet3
%     mainstream = find_mainstream(Qin,inlet,adir);
%     
%     [r1,c1]=ind2sub(size(acc), mainstream);
%     [ir,ic]=ind2sub(size(acc), inlet);
%     
%     figure
%     imagesc(log(Qin));axis image; grid on
%     hold all; plot(c1,r1,'.r','markersize',10)  %plot mainstream
%     plot(cch,rch,'.k','markersize',1)           %plot channel
%     plot(ic,ir,"o",'markersize',20)
%     title("Mainstream and Channel delineated overlaid on log(Q)")
%     legend('Mainstream - highest Q based', 'All channel', 'Mainstream Inlet')
%     c = colorbar;
%     c.Label.String = 'log(Q) (m3/s)';
% end

%% EOF
disp("###########################EOF#########################")