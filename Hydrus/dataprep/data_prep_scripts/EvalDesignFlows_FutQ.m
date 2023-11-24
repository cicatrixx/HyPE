% Evaluates design flow and related load factor for future scenarios
% Can be used for Q or Qwc
% Reads Qm%d.mat or Qm%d_wc.mat for m=1:12s
% Convert Q in m3day to m3s and get all designQ outputs in m3/s, LF: 0-1
% Saves annual avg flow (Qm13), design flow (Qdesign), annual avg load
% factor (Qdesign_LF) and mean flow through turbines given design flow (Qdesign_mean)
% Based on Qdesigner
% Script takes about 5 mins for oneQtype x 3 exceedances in HPC.

close all;
clear all;

% get paths to Q data folders
rootin = fullfile(pwd,'data','UI','data');
rootout = fullfile(pwd,'data','ASIA','Basin_UIB');

demdata="UI500m_ArcGIS.mat";
Qtype="wc"; %"nat";    %   or wc
saveQdesign=10;
% Evaluate design discharge params for following exceedances
design_exceedances=[25, 30, 40, 50, 70, 80, 90];
design_selmonths=round(design_exceedances/100*12);

addpath(genpath(fullfile(pwd,'Hydrus','dataprep')))

%% Load basin boundaries
disp('Loading basin boundaries')
load(fullfile(rootin,demdata),'outside');

%% Get input filenames and generate output folder names based on fut scenario names
fpaths_Q500m=getPaths2Files(rootin,"MLAB500m_FutQ*");  % first file listed is historical, rest are future
nQs=length(fpaths_Q500m);

for f=1:nQs
    % Get fut name from file name for all fut run
    tmps=strsplit(fpaths_Q500m{f},'MLAB500m_FutQ_m3day_');
    futname{f}=strrep(tmps{2},'.mat','');
    futout{f}=fullfile(rootout,"Fut_designQs",strcat(futname{f},'.mat'));
    % Create placeholder mat file to store Qdesigns
    fname=futname{f};
    if Qtype=="nat"
        save(futout{f},"fname")
    end
    %mkdir(fullfile(rootout,"Fut_designQs"));
end

%% Load Q_m3day file for 12 months and convert to m3/s
delete(gcp('nocreate'))
parpool(4); %not needed

parfor f=1:nQs
    Qin=fpaths_Q500m{f};
    fprintf(2,'\n\nLoading Q_%s for %s\n', Qtype,futname{f})
    Qm={};
% Load and convert Q in m3day to m3s 
    if Qtype=="nat"
        Q=load(Qin,'Q500m_m3day');
        Q500m_m3day=Q.Q500m_m3day;
        for i=1:13
            Qm{i} = Q500m_m3day{i}/(24*60*60);
        end
    elseif Qtype=="wc"
        Q=load(Qin,'Qwc_500m_m3day');
        Qwc_500m_m3day=Q.Qwc_500m_m3day;
        for i=1:13
            Qm{i} = Qwc_500m_m3day{i}/(24*60*60);
        end
    end
    [nr, nc] = size(Qm{12});


    %% Loop through all cells in matrix to get design flow vars
    for ex=1:length(design_exceedances)
        selexceedance=design_exceedances(ex); %For Q30, i.e. flow exceeded 30% of the time, David takes this as the 4th highest month, i can take percentile. QXX% is (100-XX) percentile cos exceedance is reverse of percentile
        selmon=design_selmonths(ex);
        fprintf("Processing for Qtype = %s for exceedance = %d%% \n",Qtype,selexceedance)
        % Allocate space for outputs
        Qdesign = zeros(nr,nc,'single');
        Qdesign_LF = zeros(nr,nc,'single');
        Qdesign_mean = zeros(nr,nc,'single');
        M = 4; % for running parfor in normal desktop
        for i=1:numel(Qm{12}) %
            % skip cells outside basin
            if outside(i)==1; continue; end;

            %Assigning zeros
            Qdesign2=0;
            LF_Qdesign_Monthly=0;
            LF_Qdesign=0;
            Qdesign_mean2=0;
            Qmd = zeros(12,1);

            %Collecting monthly TS data for cell
            for j=1:12; Qmd(j) = Qm{j}(i); end;

            % Qdesign is selected as XX% exceedence taken as relevant ranked month such as for 30% its the 4th highest mean monthly discharge
            Qsort = sort(Qmd,'descend');
            Qdesign2 = Qsort(selmon); %forth highest discharge month
            %Qdesign2 = prctile(Qmd,100-selexceedance);

            % Evaluate monthly and annual avg LF
            for m=1:12
                LF_Qdesign_Monthly(m) = min(1,Qsort(m)/Qdesign2)'; %LF= 1 for months where Q>=Qdesign and  LF<1 for Q<Qdesign
            end
            LF_Qdesign = mean(LF_Qdesign_Monthly);

            % Evaluate annual avg flow through the turbines
            Qdesign_mean2 = mean(min(Qdesign2,Qsort)); %Average flow rate given Qdesign

            % Archive vals to save
            Qdesign(i)=Qdesign2;
            Qdesign_LF(i)=LF_Qdesign;
            Qdesign_mean(i) = Qdesign_mean2;
        end

        %% Save all design variables for specific exceedance
        if saveQdesign
            Qdesign(outside)=0;
            Qdesign_LF(outside)=0;
            Qdesign_mean(outside)=0;

            fprintf('\nSaving vars for Qtype = %s for exceedance = %d%% \n',Qtype,selexceedance)
            % matfile_exceed = fullfile(futout{f}, sprintf('Q%d_%s.mat',selexceedance,Qtype));
            mySave1(futout{f},selexceedance,Qtype,Qdesign,Qdesign_LF,Qdesign_mean);
        end
    end

    %% Save annual avg Q
    if saveQdesign
        %matfile= fullfile(futout{f}, sprintf('Qm13_%s.mat',Qtype));
        Qm13 = Qm{13}; %Qm13
        Qm13(outside)=0;

        fprintf('Save Q_%s.mat\n',Qtype )
        mySave2(futout{f},Qtype,Qm13);
    end
end
disp('*******************************EOF*******************************')

function mySave1(filenm,Qexceedance,Qtype,Qdesign,Qdesign_LF,Qdesign_mean)
% Local save function for saving in parfor using variable names used in
% main model
%% Setup selected Qdesign Q-maps with or without water consumption
if Qtype=='wc'
    Qsuffix='_wc';
else
    Qsuffix='';
end

eval(sprintf("Q%ddesign%s = Qdesign;", Qexceedance, Qsuffix));
eval(sprintf("Q%ddesign_LF%s = Qdesign_LF;", Qexceedance, Qsuffix));
eval(sprintf("Q%ddesign_mean%s = Qdesign_mean;", Qexceedance, Qsuffix));

clearvars('Qdesign','Qdesign_LF','Qdesign_mean')
save(filenm,'-regexp',"Q..design", '-append');

end

function mySave2(filenm,Qtype,Qm13)
% Local save function for saving in parfor using variable names used in
% main model. Append to existing file
if Qtype=='wc'
    Qwc=Qm13;
    save(filenm, 'Qwc', '-append');
else
    Q=Qm13;
    save(filenm, 'Q', '-append');
end

end
