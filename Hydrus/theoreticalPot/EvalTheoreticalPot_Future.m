%% Evaluate theoretical potential for annual average dischage for all future scenarios

clear
close all
clc
addpath(genpath(fullfile(pwd,"Hydrus")))
%ncores=feature('numcores');

% Output file
save2fldr = fullfile(pwd,"output","FutRuns","fig_theory");
if ~isfolder(save2fldr); mkdir(save2fldr); end
outmat='FutTheoreticalPot_totals.mat';

% minQ threshold
minQ=0.1;%0.1m3/s
cellsz_m=500; %m
% cell spacings to evaluated theoretical potential at
dorange=ceil([500, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 25000, 40000 :20000:100000]./cellsz_m); % specified in terms of m and converted to cells...have to be multiples of 500/5000

%
modelnames=["Observed", "Historical","SSP245 GCMs", "SSP370 GCMs", "SSP585 GCMs"];
ts=[" 2036_2065", " 2066_2095"];
tname=[" Mid", " Far"];

%% Load dem based inputs
datapath=fullfile(pwd,'data','UI','data');
data='UI500m_ArcGIS.mat';
demname=fullfile(datapath,data);
fdesc=split(demname(2:end),'_');
load(demname,'fdir','dem','outside','flowdist') %,'acc','outside','catchments','basinlabels') %, 'channel') this has older channel based on facc
load(fullfile(datapath,'channel_Qbased.mat'),'channel') ;

Z_m=single(dem);
% set not basin cells to nan
Z_m(outside)=nan;

%% Get paths to future Q files
Qmon_path= fullfile(pwd,'data','UI','data');
[fpaths_Q500m, fname]=path2fldrfiles(Qmon_path,"MLAB500m*");
numScens=length(fpaths_Q500m);

% Reorder to have mid fut first and far fut later
fpaths_Q500m=fpaths_Q500m([1 2:2:25 3:2:25]);
fname=fname([1 2:2:25 3:2:25]);
fname=strrep(strrep(fname,'.mat',''),'MLAB500m_FutQ_m3day_','');

% fdetails=split(fname(2:end),'_');
%     fut_data = table(fdetails(:,1),fdetails(:,2),join(fdetails(:,4:5),"-"),repelem(["Theoretical"],24,1), ...
%         'VariableNames',{'ssp','gcm','yrs','pottype'});

%% Loop through Q files to evaluate theory pot
% Preallocate space
channel_GWh=nan(size(dorange,2),numScens);
channel_basinEnergy=[];
channel_Hgross=[];
%     channel_idxo=[];
%     channel_idx_nbr=[];
%     channel_numo=nan(size(dorange));
for f=1:length(fpaths_Q500m)
    % Get annual avg Q
    load(fpaths_Q500m{f},'Q500m_m3day','R500m_m3day');
    Q_m3s=Q500m_m3day{13}/(24*60*60);
    R_m3s=R500m_m3day{13}/(24*60*60);

    % select only Q>=1m3/s
    Q_m3sorig=Q_m3s; %select annual avg
    Q_m3s(Q_m3s<minQ)=nan;
    fprintf("Loaded Q datasets for: %s\n",fname{f})

% Archive 
Qannual_m3s_archive(:,:,f)=Q_m3s;
    %% Loop Channel Potential: Loop do and get potential
    for di=1:numel(dorange)    %For 15s res, 1 cell is ~500m. In Hydrus:50 (25km);
        % Get channel Potential for current do
        [basinEnergy, Hgross, idxo, idx_nbr]= evalChannelTheoreticalPot(Z_m, fdir, Q_m3s, dorange(di), channel, flowdist, 0);

        % Get not nan values from matrix in column vector form
        basinEnergy1 = basinEnergy(~isnan(basinEnergy(:)));
        channel_GWh(di,f)=sum(basinEnergy1(basinEnergy1(:)>=0));  %in GWh sum of not nan and not negative pot

        % Archive for saving
        if di==1
            channel_basinEnergy(:,:,f)=basinEnergy;
            channel_Hgross(:,:,f)=Hgross;
            disp("rsi=1 saved")
        end
        %channel_idxo{di}=idxo;
        %channel_idx_nbr{di}=idx_nbr;
        %channel_numo(di)=numel(idxo);
    end
    fprintf("Total channel potential at 500m in TWh in %s: %0.2f\n", fname{f},channel_GWh(1,f)/1000)
end

%% Save
% Output file
save2fldr = fullfile(pwd,"output","FutRuns","fig_theory");
if ~isfolder(save2fldr); mkdir(save2fldr); end
outmat='FutTheoreticalPot_totals.mat';

save(fullfile(save2fldr,outmat), 'channel_GWh','dorange','channel_basinEnergy','channel_Hgross','fpaths_Q500m','fname','Qannual_m3s_archive')%, ...
%'channel_idxo', 'channel_idx_nbr', 'channel_numo', 'TheorySum');

%% FINAL: Plot range of do, GWh
do_Gernaat=25; % in km
do_Hoes= 225/1000; % in km (pg 3)
modelnames=["Observed", "Historical","SSP245 GCMs", "SSP370 GCMs", "SSP585 GCMs"];

figure
c_ssp=[mycbrewer('Greens',4);mycbrewer('Blues',4);mycbrewer('YlOrBr',4)];
colororder(c_ssp)
hold all
l1=plot(dorange*cellsz_m/1000, channel_GWh(:,2:11)/1000,"-","LineWidth",1); %only mid fut,'DisplayName','Foresight-based')%'Color',ccol(1,:),
l2=plot(dorange*cellsz_m/1000, channel_GWh(:,12:24)/1000,"--","LineWidth",1);%only far fut,'DisplayName','Foresight-based')%'Color',ccol(1,:),

l3=plot(dorange*cellsz_m/1000, channel_GWh(:,1)/1000,".-k","LineWidth",2); %historical ,'DisplayName','Hindsight-based')%'Color',ccol(1,:),
ylabel('Total energy (TWh/yr)')
xlabel('River spacing (km)')
title(sprintf("Total theoretical potential for Q>=%0.1fm^3/s",minQ))
grid on
legend([l1([4:4:12]); l2([4:4:12]); l3], [strcat(modelnames(3:end), " mid-future") strcat(modelnames(3:end), " far-future") "Historical"])


%orient('landscape')
%print('G:/GitConnect/output/TheoryPot_WrapUp/TheoryPot_wSpacing.pdf','-dpdf','-fillpage')

