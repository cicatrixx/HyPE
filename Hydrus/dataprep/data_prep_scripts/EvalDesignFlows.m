% Evaluates design flow and related load factor. Can be used for Q or Qwc
% Reads Qm%d.mat or Qm%d_wc.mat for m=1:12s
% All Q outputs in m3/s, LF: 0-1
% Saves annual avg flow (Qm13), design flow (Qdesign), annual avg load
% factor (Qdesign_LF) and mean flow through turbines given design flow (Qdesign_mean)
% Based on Qdesigner
% Script takes about 5 mins for oneQtype x 3 exceedances in HPC. 

close all;
clear all;
root_data = fullfile(pwd,'data');
% get paths to continent data folders
root = fullfile(root_data,'UI','data');

Qin="MLAB500m_40yrClimatology_m3day_R4.mat";
demdata="UI500m_ArcGIS.mat";
Qtype="wc"; % "nat";%  or wc

% Evaluate design discharge params for following exceedances
design_exceedances=[25, 30, 40, 50, 70, 80, 90];
design_selmonths=round(design_exceedances/100*12);

%% Load Q file for 12 months and convert to m3/s
disp('Loading Q')
if Qtype=="nat"
    load(fullfile(root,Qin),'Q500m_m3day');
    for i=1:13
        Qm{i} = Q500m_m3day{i}/(24*60*60);
    end
elseif Qtype=="wc"
    load(fullfile(root,Qin),'Qwc_500m_m3day');
    for i=1:13
        Qm{i} = Qwc_500m_m3day{i}/(24*60*60);
    end
end
[nr, nc] = size(Qm{12});

%% Load basin boundaries
disp('Loading basin boundaries')
load(fullfile(root,demdata),'basinmask','catchments','outside');

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
    fprintf('\nSaving vars for Qtype = %s for exceedance = %d%% \n',Qtype,selexceedance)
    matfile_exceed = fullfile(root, sprintf('Q%d_%s.mat',selexceedance,Qtype));
    save(matfile_exceed,'Qdesign','Qdesign_LF','Qdesign_mean');
end

%% Save annual avg Q
matfile= fullfile(root, sprintf('Qm13_%s.mat',Qtype));
Qm13 = Qm{13}; %Qm13
disp(sprintf('Save Q_%s.mat',Qtype ))
save(matfile,'Qm13');

disp('*******************************EOF*******************************')
% end

%% Check
% figure(1);clf;imagesc(Qdesign_LF);axis image
% r = 817;
% c = 2256;
%
% %Collecting monthly data
% for j=1:12; Qmd(j) = Qm{j}(r,c); end;
%
% % Qdesign is allows 30% exceedence, so 70% of maximum flow
% % Here we take the 4th highest mean monthly discharge
% Qsort = sort(Qmd,'descend');
% Qdesign2 = Qsort(4) %forth highest discharge month

