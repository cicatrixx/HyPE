%% Route 5km R time series .tiff files to get discharge for future scenarios
% Created By    : Sanita Dhaubanjar
% Created For	: SustaIndus project WP2
% uses parfor!
% For time series, we use the 
% R in mm/MONTH to generate Q and Qwc in m3/MONTH
% SHOULD UPDATE CODE NEXT TIME TO SAVE M3/DAY AND ONLY SAVE BASIN BOUNDS DIRECTLY - 23 Dec 2023

%hist run takes 48GB RAM
clc
clear all
close all

nworkers =4; %cant do 4 in 64gbHPC                     % M specifies maximum number of parpool workers
root=pwd;

matpath=fullfile(root,"data","UI","data","Q_monTS_mmmonth");
saveQ2mat=fullfile(matpath,"MLAB500m_FutQ_TS_m3day_");  % add futname as suffix here
createplots=1;
nodata32bitFlo = -340282346638528859811704183484516925440.00; % for 32 bit signed FLOAT type
showplot=0;
setnan=0;

%input paths
root5km_TS_tiff=fullfile(pwd,"data","data_prep","UIB_outputs_monTS_mmmonth");
addpath(genpath(fullfile(pwd,'Hydrus')))

%% Load 500m basin data acc, adir
disp('Loading acc and dir at 500m');
data500m='UI500m_ArcGIS.mat';
load(fullfile(matpath,data500m), 'acc', 'adir')%,'outside','basinmask','channel') %,  'fdir', 'outlet'

%% Create flist
tifflist=compose("TotrM%05d.tiff",1:12*40)';
fpaths=path2fldrfiles(root5km_TS_tiff,"monTS*");
fdetails=split(fpaths(2:end),filesep); % Just take future file details. idx 1 is hist
n=size(fdetails,2);
fname=["hist"; join(fdetails(:,(n-2):n),'_')];

disp(['Created list of folders w tiff files']);

%% Get Q natural
    delete(gcp('nocreate'));
    parpool('local',nworkers);
parfor (f= 1:9,nworkers)
    %% Make Rm by loading .tiff files
    Rm_SPHY_mmmonth={};
    if f==1 %for historical runs %40 yrs data
        seltiff=12*40; %40 yrs
        
    else  % for future runs
        seltiff=12*30;
    end
    for m=1:length(tifflist(1:seltiff))
        Rm_SPHY_mmmonth{m} = loadSPHYtiff(fullfile(fpaths{f},tifflist{m}), nodata32bitFlo , showplot,setnan);
    end

    disp(['Loaded TotR tiff  files for scenario: ' fpaths{f}]);
        
    %% get Q natural
    fut_name=fname{f};
    scale_mm_to_m3day = 500^2*1e-3;  %for 500m cell size
% THOUGH THIS SAVES IN M3DAY, the unit is actually M3/MONTH! i did not rerun script again so script is not changed here
    [Q500m_m3day, R500m_m3day]=downscaleR2Q(Rm_SPHY_mmmonth, acc, adir, scale_mm_to_m3day);

    fprintf('Routed Qnatural for %s\n',fut_name);




    %% Save mat file
    mySave(strcat(saveQ2mat,fut_name,'.mat'),Q500m_m3day, R500m_m3day)%, 'Qwc_500m_m3day', 'Rwc_500m_m3day');
    disp('Qnatural and Qwc files saved!');
end

%% EOF
disp("###########################EOF#########################")

%%
function mySave(filenm,Q500m_m3day, R500m_m3day)
% Local save function for saving in parfor using variable names used in
% main model. Append to existing file

    save(filenm, '-v7.3','Q500m_m3day', 'R500m_m3day');

end

