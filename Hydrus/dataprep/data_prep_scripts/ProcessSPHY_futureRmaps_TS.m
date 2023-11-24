%% Save SPHY and PCRASTER based .map files for monthly time series TotR as .tiffs
% Generate input and output filenames for user specified timeframes.
% Loads .maps and cover to .tiffs
% Loop through filepaths for all future scenario folders inside root5km
% All R unit is mm/MONTH because we use monthly flux files!
% USES GDAL_TRANSLATE for map2tiff coversion!
% Created By    : Sanita Dhaubanjar on 16 Nov 2022
% Created For	: SustaIndus project WP2
%=========================
clc
close all
clear all

%output paths
root5km_TS_tiff=fullfile(pwd,"data","data_prep","UIB_outputs_monTS_mmmonth");

%input paths
root5km_TS =fullfile(strrep(pwd,"HPmodel",""),"SPHY_UIB","UIB_outputs_all");
addpath(genpath(fullfile(pwd,'Hydrus\')))

%% Get paths to all folders insider root
fpaths_in=path2fldrfiles(root5km_TS);

%% Generate file indices and filenames
datemin = [1979 2015 2015];
datemax = [2018 2100 2100];
startyr_lst= [1979 2036, 2066];
endyr_lst= [2018 2065, 2095];

for tf=1:3
    nyears=endyr_lst(tf)-startyr_lst(tf)+1;
    nmonths = 12*nyears;

    % Generate monthly time series from Feb 1 to Jan 1
    mon_dates = datetime(startyr_lst(tf),1,1) + calmonths(1:nyears*12)';
    mon_fidx= datenum(mon_dates) - datenum(datemin(tf),1,1);

    mon_fname_in{tf}=compose("TotrM%07.3f",mon_fidx/1000);
    mon_fname_out{tf}=compose("TotrM%05d.tiff",1:length(mon_fname_in{tf}));
end
disp("Created input output filenames")

%% Loop through historical scenario
tf=1;
f=1;
disp(['Loading TotR from map files for scenario: ' fpaths_in{f}]);

parfor m = 1:length(mon_fname_in{tf})
    %% Convert TotR .map  files to .tif
    fpaths_out=fullfile(strrep(fpaths_in{f},root5km_TS, root5km_TS_tiff),sprintf("monTS_%d-%d",startyr_lst(tf),endyr_lst(tf)));
    if m==1; mkdir(fpaths_out); end
    R_5kmtif=convertmap2tif(fullfile(fpaths_in{f},mon_fname_in{tf}(m)),  fullfile(fpaths_out,mon_fname_out{tf}(m)));
end

%% Loop through filepaths for all future scenario, load 5km runoff maps or open preloaded one
for f = 2:length(fpaths_in)
    sspfpaths=path2fldrfiles(fpaths_in{f});
    for cpath_in=sspfpaths'
        disp(['Loading TotR from map files for scenario: ' cpath_in{:}]);

        for tf=2:3
            cpath_out=fullfile(strrep(cpath_in{:},root5km_TS, root5km_TS_tiff),sprintf("monTS_%d-%d",startyr_lst(tf),endyr_lst(tf)));
            
            parfor m = 1:length(mon_fname_in{tf})
                if m==1; mkdir(cpath_out); end

                %% Convert TotR .map  files to .tiff
                R_5kmtif=convertmap2tif(fullfile(cpath_in{:},mon_fname_in{tf}(m)),  fullfile(cpath_out{:},mon_fname_out{tf}(m)));
            end
        end
    end
end

%% EOF
disp("###########################EOF#########################")

