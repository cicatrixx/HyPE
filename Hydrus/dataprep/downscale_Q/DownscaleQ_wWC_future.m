%% Route preloaded 5km R to get discharge w and wo water consumption for future scenarios
% downscaleQ_wWC.m
% Created By    : Sanita Dhaubanjar on 24 Nov 2020
% Created For	: SustaIndus project WP2
% uses parfor!
% in mm/day to generate Q and Qwc in m3/day

clc
clear all
close all

nworkers =0;                     % M specifies maximum number of parpool workers
root=pwd;

matpath=fullfile(root,"data","UI","data");
saveQ2mat=fullfile(matpath,"MLAB500m_FutQ_m3day_");  % add futname as suffix here
createplots=1;
addpath(genpath(fullfile(pwd,'Hydrus\')))

%% Load 500m basin data acc, adir
disp('Loading acc and dir at 500m');
data500m='UI500m_ArcGIS.mat';
load(fullfile(matpath,data500m), 'acc', 'adir','outside','basinmask','channel') %,  'fdir', 'outlet'

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

%% Load 5km runoff maps or open preloaded one
saveR2mat='SPHY5km_FutureQ_mmday.mat';
load(fullfile(matpath,saveR2mat),'Rm'); % Rm has data for all 24 scenarios.

%% get Q natural
nfuts=length(Rm.name);
for (f= 1:nfuts)
    fut_name=Rm.name{f};
    Rm_SPHY_mmday=Rm.data{f};

    %% get Q natural
    scale_mm_to_m3day = 500^2*1e-3;  %for 500m cell size
    [Q500m_m3day, R500m_m3day]=downscaleR2Q(Rm_SPHY_mmday, acc, adir, scale_mm_to_m3day);
    fprintf('Routed Qnatural for %s\n',fut_name);
%     % Archive all 500m Q and R
%     R500m{f}=R500m_m3day;
%     Q500m{f}=Q500m_m3day;
%     % Saved annual Q and R to geotiff
    savemat2Pantpetiff(fullfile(root,"data","data_prep","UIB_outputs_LTMavgs",strcat(fut_name,"_Q500m_M13_m3s.tif")), Q500m_m3day{13}/(24*60*60))
    savemat2Pantpetiff(fullfile(root,"data","data_prep","UIB_outputs_LTMavgs",strcat(fut_name,"_R500m_M13_m3s.tif")), R500m_m3day{13}/(24*60*60))

    %% get Qwc
    [Qwc_500m_m3day, Rwc_500m_m3day]=downscaleR2Q(Rm_SPHY_mmday, acc, adir, scale_mm_to_m3day,WCt_mmday);
    fprintf('Routed Qwc for %s\n',fut_name);
    % Archive all 500m Q and R
%     Rwc_500m.data{f}=Rwc_500m_m3day;
%     Qwc_500m.data{f}=Qwc_500m_m3day;

    %% Save mat file
    save(strcat(saveQ2mat,fut_name,'.mat'),'-v7.3','Q500m_m3day', 'R500m_m3day', 'Qwc_500m_m3day', 'Rwc_500m_m3day');
    disp('Qnatural and Qwc files saved!');
end


%% Save mat file
%save(saveQ2mat,'-v7.3','Q500m', 'R500m');%, 'Qwc_500m', 'Rwc_500m');
%disp('Qnatural and Qwc files saved!');

%% EOF
disp("###########################EOF#########################")