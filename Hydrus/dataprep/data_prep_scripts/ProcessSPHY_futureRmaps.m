%% Compile SPHY and PCRASTER based .map files for monthly long term TotR into one .mat file for future scenarios
% Read 5km .map SPHY outputfiles for TotR (for 12 months). Loads .maps as .tiffs then save into .mat files
% Generate annual TotR as M13 by average monthly files.
% Loop through filepaths for all future scenario folders inside root5km
% All R unit is mm/day
% USES GDAL_TRANSLATE for map2tiff coversion!
% Created By    : Sanita Dhaubanjar on 16 Nov 2022
% Created For	: SustaIndus project WP2
%=========================
clc
close all
clear all

%output paths
matpath=fullfile(pwd,'data/UI/data');
saveR2mat='SPHY5km_FutureQ_mmday.mat';

%input paths
%root5km =fullfile(pwd,"..","\SPHY_UIB\UIB_outputs_LTMavgsssp245\NorESM2-MM\LTavgs_2036_2038");
root5km =fullfile(pwd,"data\data_prep\UIB_outputs_LTMavgs");
addpath(genpath(fullfile(pwd,'Hydrus\')))
sspls= {'ssp245', 'ssp370', 'ssp585'}; % dir(root5km);

check_annualsum=1;
%% Get paths to all LTavgs folders
tmp=dir(fullfile(root5km, '**/LTavgs*'));
fpaths=join([{tmp.folder}', {tmp.name}'],filesep,2);
fprintf("\n\n\n\nLoading map files for %d scenarios\n", length(fpaths) )

%% Prep filenames for 12 months
for m=1:9
    monfname(m)=sprintf("TotrA%dM0.360",m);
end
for m=10:12
    monfname(m)=sprintf("TotrA%dM.360",m);
end

%% Loop through filepaths for all future scenario, load 5km runoff maps or open preloaded one
nodata32bitFlo = -340282346638528859811704183484516925440.00; % for 32 bit signed FLOAT type
showplot=0;
setnan=0;
fidx=1;
for cpath = fpaths'
    % Get current future scenario name from current filepath
    ctmp=strsplit(cpath{:},'\');
    scname=strjoin(ctmp(7:9),'_');

    disp(['Loading TotR from map files for scenario: ' scname]);
    Rm_SPHY_mmday={};
    for m = 1:12
        %% Convert TotR .map  files to .tiff
        R_5kmtif=convertmap2tif(fullfile(cpath{:},monfname{m}),  fullfile(cpath{:},strcat('TotrA',num2str(m),".tiff")));
        Rm_SPHY_mmday{m} = loadSPHYtiff(R_5kmtif, nodata32bitFlo , showplot,setnan);
    end
    % Prepare annual Totr
    Rm_SPHY_mmday{13} =  sum(cat(3,Rm_SPHY_mmday{1:12}),3)/12;

    % Archive all future Rs
    Rm.name{fidx}=scname;
    Rm.data{fidx}=Rm_SPHY_mmday;
    fidx=fidx+1;
end


%% Save to .mat
disp('Saving Rm data')
matfile=fullfile(matpath,saveR2mat);
save(matfile,'-v7.3','Rm');

%% Plot historical and future R and apply labels and formatting to graph
if check_annualsum
    % Plot monthly sum of R over the whole basin
    figure(10);clf; hold on
    colororder(cbrewer2('Paired',8))
    %colororder([cbrewer2('Greens',8); cbrewer2('Purples',8);cbrewer2('Oranges',8) ]) % 2LTs x 3ssp x 4models
    
    for f=1:length(fpaths)
        Rm_SPHY_mmday=Rm.data{f};
        % Plot monthly sum of R over the whole basin
        plot(squeeze(sum(sum(cat(3,Rm_SPHY_mmday{:})))),'DisplayName',Rm.name{f})
    end

    % Add historical R
    load(fullfile(matpath,'SPHY5km_40yrClimatology.mat'),'Rm_SPHY_mmday');
    plot(squeeze(sum(sum(cat(3,Rm_SPHY_mmday{:})))),'o-k','DisplayName',"Historical", "LineWidth",1.5)

    % Apply labels and formatting to graph
    legend('Interpreter','none')
    box on
    title("Spatial sum of 5km monthly average runoff maps for different future scenarios")
    set(gca,'xtick',1:13,'linestyleorder', {':o','-.s','-d'},...
        'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','ANNUAL'})
    xline(12.5,'LineStyle',':')
    ylabel("SPHY 5km Runoff in mm/day")  

end
