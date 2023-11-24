addpath(genpath(fullfile(pwd,"..","Hydrus")))

% varname
nbasin=103;
continent_in="ASIA";
runnameidx=103; %num in filenames
r_prefix=sprintf('FutR%d_',runnameidx);

rootf=pwd; %path to HPmodel 'G:\SurfDrive\HPmodel'; %
rootif=fullfile(rootf,sprintf('\\data\\%s\\Basin_UIB\\PantpeBasin_%d.mat',continent_in,nbasin));
rootof=fullfile(rootf,'output');
rootoffut=fullfile(rootf,'output','FutRuns_Figs_Analysis');


%% params
costlim=0.1; %$0.01/kWh COE limit COEAlls<=costlim
cellsz_m=500; %m
energyreq_MWhpercapita=0.6; %MWh/capita --not convert to TWh so out pop is in millions

%% Read model order
modorder=readtable(fullfile(rootof,"FutRuns_Figs_Analysis" ,"FutScen_Tracker.xlsx"),Sheet="GCMdetails", Range="A1:F13");

%% stringnames
tframe=["2036-2065", "2066-2095"];
tname=["Mid","Far"];
termname=["Medium","Long"];
%
sspshort=[ "ssp245" "ssp370" "ssp585"];
%sspnames=["SSP 2 + RCP 4.5" "SSP 3 + RCP 7.0" "SSP 5 + RCP 8.5"]; %AL proposed
sspnames=["RCP 4.5 + SSP 1" "RCP 7.0 + SSP 2" "RCP 8.5 + SSP 3"]; %Mine
rcpnames=["RCP 4.5" "RCP 7.0" "RCP 8.5"];
tf_full=strcat(tname,": ", tframe);
cornernames=["Warm Wet", "Warm Dry","Cold Dry", "Cold Wet"];
cornernamesshort=["WW","WD","CD","CW"];
modelnames=["Observed", "Historical",rcpnames];
%
pottypes5={'Theoretical','Technical','Financial','Sustainable','Visualized'};
pottypes5_short={'Theory','Tech','Fin','Sust','Vis'};
pottypes6=[pottypes5(:)', {'Existing'}];
pottypes3={'Technical','Financial','Sustainable'};
pottypes3_short={'Tech','Fin','Sust'};
pottypes2=pottypes3([1,3]);
% w sust-tech
newpots4=["Tech-Tech" "Tech-Fin" "Sust-Tech" "Sust-Fin"];
newpots4cl=["Technical","Financial","Sust-Tech","Sustainable"];
%
planttypes=["River Power Plant","Diversion Canal Plant"];

%% rcp ssp cornernames
rcpcornernames=[];
ssprcpcornernames=[];
rcp1cornernames=[];
cornernamesrcp1=[];
for i=1:3
    rcp1cornernames=[rcp1cornernames strcat(rcpnames(i),": ",cornernames(1)) cornernames(2:end)];
    rcpcornernames=[rcpcornernames strcat(rcpnames(i),": ",cornernames)];
    ssprcpcornernames=[ssprcpcornernames strcat(sspnames(i),": ",cornernames)];
    cornernamesrcp1= [cornernamesrcp1 strcat(rcpnames(i),": ",cornernames(1)) cornernames(2:end)];
end
histrcpcornernames=["Historical" rcpcornernames rcpcornernames];
histshortcornernames=["Historical" cornernamesshort cornernamesshort cornernamesshort];
cornernamesrep=repmat(cornernames,1,3);
% add TF to names
tfrcpnames=[strcat(tname(1)," ",rcpnames(1)) rcpnames(2:3) strcat(tname(2)," ",rcpnames(1)) rcpnames(2:3)];
tfhistrcpnames=["Historical: 1979-2018" tfrcpnames];
tfhistrcpcornernames=["Historical: 1979-2018" strcat(tname(1)," ",rcpcornernames(1)) rcpcornernames(2:end) strcat(tname(2)," ",rcpcornernames(1)) rcpcornernames(2:end)];

%% Timeframes
hist_tframe='1979-2018';
hist_yrs=40;
fut_yrs=30;

%annual TS
histyrs_ts=1979:2018;
fyrs_ts1=2036:2065;
fyrs_ts2=2066:2095;

%monthly TS
%histmon_ts = datetime('January 1, 1979') + calmonths(0:12*hist_yrs-1)';
histmon_ts = datetime(histyrs_ts(1),1,1) + calmonths(0:12*hist_yrs-1)';
futmon_ts1 = datetime(fyrs_ts1(1),1,1) + calmonths(0:12*fut_yrs-1)';
futmon_ts2 = datetime(fyrs_ts2(1),1,1) + calmonths(0:12*fut_yrs-1)';

%% Subbasin labels
load(fullfile(rootf,'data\UI\data\UI500m_ArcGIS.mat'), 'basinlabels','outside')
basindata=basinlabels(1:8,1:3);
label_subbasin_basin =[basindata.basinnames{:}, "All of upper Indus"];

%% Fig things
baralpha=0.7;
energylabel='Energy (TWh yr^{-1})';%'Energy in TWh/yr';
%text(xx(1:4),repelem(11800,4),tot_allscen.pottype(1:4:16),'HorizontalAlignment','center','FontAngle','italic')
  %  ylabel('MWh per capita per year','fontweight','bold')
charlbl =  compose("(%s) ",('a':'i').'); % {'(a)','(b)','(c)','(d)'}


%% Color palette for fut scenarios
c_ssp=[mycbrewer('Greens',4);mycbrewer('RdPu',4);mycbrewer('YlOrBr',4)]; % colors for 3 ssps x 4 corners
c_ssp_main=[brighten(c_ssp(4:4:12,:),0.5); brighten(c_ssp(4:4:12,:),0.5)];
%c_ssp_main=cbrewer2('Dark2',4);

c_ssp26=[0 0 0; c_ssp; c_ssp; 1 1 1]; % add black for historical, add white for white space
%% Other color for plots
cl_RP= [216,179,101]/255; %brownish %[112 173 71]/255; %'#70AD47';
cl_DP= 127/255*[1 1 1]; %gray %'#7F7F7F';
cl_channel= [.53 .81 1]; % very light river %[0 0.45 0.74]; %river blue

mycolors={cl_RP cl_DP};
baralpha=0.7;
myalpha=0.6; % for scatter
mygraylines=.15*[1 1 1]; % for axes borders and boxes
mybrown=[165, 20, 20]/255;
subbasinalpha=0.9;
extraaxis=[1 0.1 0.2]; % light red
%% subbasin colors: from Wong
%from: https://www.nature.com/articles/nmeth.1618/figures/2
cmap8_wong_gray= [
    230, 159, 0 %Orange
    86,180,233 %Sky Blue
    0, 158, 115 %Bluish green
    240, 228, 66 %Yellow
    0, 114, 178 %Blue
    213, 94, 0 %Vermillion
    200, 200, 200 %Grey
    204,121,167 %Reddish purple
    ]/255;

%% future potential colors
c_pot5=cbrewer2('Set2',5);
% reorder for sust to be green
c_pot5=c_pot5([4,2,3,1,5],:);
%% For each ssp and Q, get idx for all pfs related to it (4 corners for each ssp idx)
startidx=2:4:25; % for ssps
endidx=5:4:25;
% for histQ
selpf4qf{1}=1:25;

for i=1:3*2 % 3sspx2tf
    % Get corners for each ssp
    sspgrps{i}=startidx(i):endidx(i); % contains 3ssps for near fut followed by far fut
    % Repeat sspgroup for each corner Q in that sspgroup
    for ii=sspgrps{i}
        selpf4qf{ii}=sspgrps{i};
    end
end

%% Figure settings
subplottight = @(m,n,p) subtightplot (m, n, p, [0.005 0.005],[.08 .08],[.08 .01]); %[vgap hgap], marg_h -[low up],marg_w -[l r]

%%
%xline(51,'k','LineWidth',2) % for tf
%xline([17 34 68 85],'Color',mygraylines,'LineWidth',1.15) % for ssp