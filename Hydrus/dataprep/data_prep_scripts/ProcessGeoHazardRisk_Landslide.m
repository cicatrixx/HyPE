% ProcessGeoHazardRisk.m
% Created By    : Sanita Dhaubanjar on 27 Oct 2021
% Created For	: SustaIndus project WP2
%=========================
% Code for processing Landslide Susceptibility files

close all
clear
root_data= fullfile(pwd,'\data\data_prep');
suffix="_UIprj"; % suffix added to output filename
matpath=fullfile(pwd,'\data\UI\data');
omat='LandSlideSusceptibility_R2.mat';
addpath(genpath(fullfile(pwd,'\Hydrus\devFiles'))) % add all subfolders within folder

genLocalClass=0;

%% Load .mat with LS fuzzy values 0-1 in PCS
load(fullfile(matpath,omat),'data');
load(fullfile(pwd,'data\UI\data\UI500m_ArcGIS.mat'), 'outside')

%% Crop to basin and sort
LS_fuzzy=maskBasin(data,~outside);

% Sort from high to low vals - only vals in basin to avoid nans
LS_fuzzy_sort = sort(LS_fuzzy(outside==0),'descend');

%% Evaluate local risk classes or load global ones
if genLocalClass
    %% Find index number for breakpoints
    outsuffix='UIB_Local_classified';
    nclass =5;
    r=2;
    nbins=r.^(0:nclass-1); %class definition is based on a geometric progression whereby consecutive class has 2(=r)x the # of cells in current class
    pop=numel(LS_fuzzy_sort);
    pop_perbin=nbins*round(pop/sum(nbins));
    for i=1:nclass-1
        breakpoint_idx(i)=sum(pop_perbin(1:i));
    end
    breakpoint_idx(nclass)=numel(LS_fuzzy_sort); % first boundary should be lowest val in data but because of roundind errors some cells are being missed out
    breakpoints = LS_fuzzy_sort(breakpoint_idx);
    breakpoints =flipud(breakpoints); % Reorder from VLow to Vhigh class
else
    % Global risk classes from Thomas's email
    outsuffix='UIB_Global_classified';
    breakpoints =[0 0.1038620 0.4383037 0.6125923 0.6864969];
end

%% Plot sorted data and breakpoints
figure;
plot(LS_fuzzy_sort)
hold on
xylines(breakpoints,':r','Y');

%% Assign hazard classes based on local or global risk classes
landslide_level= zeros(size(LS_fuzzy));
lbounds= reshape(breakpoints,[numel(breakpoints),1]); % read breakpoints as column vector
ubounds= [lbounds(2:end); 1];

for i=1:length(lbounds)
    landslide_level(LS_fuzzy>=lbounds(i)& LS_fuzzy<ubounds(i))=i; 
end
%%
figure;
subplot(1,2,1)
imagescnan(landslide_level)
title(outsuffix)
subplot(1,2,2)
histogram(landslide_level(landslide_level>0))
grid on
xticklabels({'Very low','Low','Med','High','Very high'})


%% Save LS_hazard to .tif
save(fullfile(matpath,['LandSlideLevels_',outsuffix,'.mat']),'landslide_level')
landslide_level2=maskBasin(landslide_level,~outside,-999);
savemat2Pantpetiff(fullfile(root_data,'LHASA_UpdateFromStanley2021_CentralAsia',[outsuffix,'_UIprj.tif']),landslide_level2)
disp("***************************************EOF***************************************")

