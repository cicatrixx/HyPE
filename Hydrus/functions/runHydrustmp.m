function runname= runHydrustmp(runname_prefix, rootpath,...
    scenario, policytype, runsustain, geohazard_scentype, geohazard_select,...
    makechanges2main, makechanges2cost)
% full code setup for running the UIB hydropower potential exploration model
% Mandatory Inputs:
    % runname_prefix- Run name is created by adding this prefix + suffixes indicating Full/Rem,Large/med/mix, TechFin/Sust runs and geohazard scen
    % rootpath      - main directory with the Hydrus, data and output folders
    % scenario      - 'Full' or 'Remain'
    % policytype    -  % policy scenarios 0=trial w mixed, 1= trial w large, 2= large focus, 3= medium focus, 4= mixed focus
    % runsustain    - (0/1) turn off/on sustainable constrains;
    % geohazard_scentype - % For sustainable runs, 0=no geohazard consideration, 1=Status quo, 2= Risk averse, 3= Multi hazard
    % geohazard_select-[1 1 1 1]; %geohazard data to use Apply EQ-PGA, EQ-Thrust, LS, GLOF hazard risks
% Optional inputs:
    % makechanges2main - string to be fed into eval function to make changes to vars defined in main function
    % makechanges2cost - string to be fed into eval function to make changes to vars defined in costconfig
   

%RP=River power Plants and DP=Diversion canal plants
tic
close all

% Renamed split_main_trib to n_river_level
continent_in = 'ASIA';
CID= 2; % previously this was inside georefs
nbasin = 101; % basin file number to load
savedata=1;
showcheckplots=0;

% vars to change for sensitivity analysis: 
%Main function: sel_design_flows, storage_lim, Plarge_minH, Psmall_minH=20, slackflow
%Cost config: eta_generation_DP_large, eta_generation_DP_small, eta_generation_RD, interest
%% String parameters
policynames={'TrialwMix','TrialwLarge','Large','Medium','Mixed'};
geohazard_scennames={'NoGeohazard', 'CostBased','RiskAverse', 'Multi-hazard'};
geohazard_snames={'EQ','EQ', 'LS', 'GLOF'};  % short ones used to gen filename

%% Setup policy scenario settings
if policytype==0 %TRIAL RUN W MIXED
    do_main_km=1000;            %For main Distance between the outlets (spacing), in km Standard:(4km) set to nan if not used
    do_trib_km=500;            %For tributary Distance between the outlets (spacing), in km Standard:(2km); for development: 2000(1000km)
    minDis_main_km=990;           %Minimum distance between cascading plants
    minDis_trib_km=498.5;         %Minimum distance between cascading plants
    varyRiverLevel_HPtype=2;              %Use different plant type and river spacing depening on stream type: (0=)use same spacing and RP/DP in all (1=) mainstream (RP/DP) and tributaries (only DP) (=2) primary (RP/DP), secondary (DP), tertiary (DP)
    enablesmallPdesignQH = 1;       %Use different Qdesign and add smallP_minH threshold
    enablesmallPcost = 1;           %Use different cost function and add minH threshold
    
    %savedata=0;
elseif policytype ==1 %TRIAL RUN W LARGE
    do_main_km=1000;            %For main Distance between the outlets (spacing), in km Standard:(4km) set to nan if not used
    do_trib_km=do_main_km;      %River segment length for tributary Distance between the outlets (spacing), in km Standard:(2km); for development: 2000(1000km)
    minDis_main_km=990;           %Minimum distance between cascading plants
    minDis_trib_km=minDis_main_km;         %Minimum distance between cascading plants
    %savedata=0;
    
    varyRiverLevel_HPtype=0;              %Use different plant type and river spacing depening on stream type: (0=)use same spacing and RP/DP in all (1=) mainstream (RP/DP) and tributaries (only DP) (=2) primary (RP/DP), secondary (DP), tertiary (DP)
    enablesmallPdesignQH = 0;       %Use different Qdesign and add smallP_minH threshold
    enablesmallPcost = 0;           %Use different cost function and add minH threshold
    
elseif policytype==2 % LARGE FOCUS
    do_main_km=25;               %For main Distance between the outlets (spacing), in km Standard:(4km) set to nan if not used25
    do_trib_km=do_main_km;      %River segment length for tributary Distance between the outlets (spacing), in km Standard:(2km); for development: 2000(1000km)
    minDis_main_km=1;           %Minimum distance between cascading plants
    minDis_trib_km=minDis_main_km;         %Minimum distance between cascading plants
    varyRiverLevel_HPtype=0;              %Use different plant type and river spacing depening on stream type: (0=)use same spacing and RP/DP in all (1=) mainstream (RP/DP) and tributaries (only DP) (=2) primary (RP/DP), secondary (DP), tertiary (DP)
    
    enablesmallPdesignQH = 0;       %Use different Qdesign and add smallP_minH threshold
    enablesmallPcost = 0;           %Use different cost function and add minH threshold
elseif policytype==3 % MEDIUM FOCUS
    do_main_km=4;               %River segment length for main Distance between the outlets (spacing), in km Standard:(4km) set to nan if not used
    do_trib_km=do_main_km;      %River segment length for tributary Distance between the outlets (spacing), in km Standard:(2km); for development: 2000(1000km)
    minDis_main_km=1;           %Minimum distance between cascading plants
    minDis_trib_km=minDis_main_km;         %Minimum distance between cascading plants
    varyRiverLevel_HPtype=0;              %Use different plant type and river spacing depening on stream type: (0=)use same spacing and RP/DP in all (1=) mainstream (RP/DP) and tributaries (only DP) (=2) primary (RP/DP), secondary (DP), tertiary (DP)
    
    enablesmallPdesignQH = 0;       %Use different Qdesign and add smallP_minH threshold
    enablesmallPcost = 0;           %Use different cost function and add minH threshold
elseif policytype==4 % MIXED FOCUS
    do_main_km=4;               %River segment length for main Distance between the outlets (spacing), in km Standard:(4km) set to nan if not used
    do_trib_km=2;               %River segment length for tributary Distance between the outlets (spacing), in km Standard:(2km); for development: 2000(1000km)
    minDis_main_km=1;           %Minimum distance between cascading plants
    minDis_trib_km=0.5;         %Minimum distance between cascading plants
    
    varyRiverLevel_HPtype=2;              %Use different plant type and river spacing depening on stream type: (0=)use same spacing and RP/DP in all (1=) mainstream (RP/DP) and tributaries (only DP) (=2) primary (RP/DP), secondary (DP), tertiary (DP)
    enablesmallPdesignQH = 1;       %Use different Qdesign and add smallP_minH threshold
    enablesmallPcost = 1;           %Use different cost function and add minH threshold
end

%% Sustain constraints - yes(1) or no(0) - the ones that may be changed
if runsustain
    waterconsumption    =1; %Apply water consumption
    slackflow_constraint=1; %Enable slackflow of x% to ensure natural flow of the river (see slackflow -setting)
    protarea_constr     =1; %Apply no dams in protected areas
else
    waterconsumption    =0; slackflow_constraint=0;
    protarea_constr     =0;
    geohazard_constr    =0; geohazard_cost=0;
end

%% Settings that are changed
slackflow=30/100;           %Percentage of water diverted around the dam to ensure natural ecological flow of the river
enablereslimit = 1;         %Restrict reservoir storage volume to be a % of annual average inflow volume
storage_lim=5/100;          %Maximum limit on reservoir size as percentage of annual average inflow volume
Plarge_minH=4;              %Minimum head in m for large P (default =4m)
Psmall_minH=20;             %Minimum head in m for small P (default =20m)
RD_minH=4;                  %Minimum dam height or head for RD project in m
sel_design_flows = {'Q40', 'Q30', 'Q80'};   %design flow levels for RP, large DP and small DP

%% Make changes to config vals if any
if exist('makechanges2main','var')
    disp("Changes made to mainconfig")
    eval(makechanges2main);
else
    makechanges2cost='';
end

%% Setup Full/Remain scenario
if strcmp(scenario,'Full')
    ExistingDams_constraint=0;
elseif strcmp(scenario,'Remain')
    ExistingDams_constraint=1;%No dams on locations in existing hydro lakes
end

%% Fixed Settings
% RD = River power systems (River dams)
% DP = Reservoir systems (Damp-Pipe systems)
% P  = Diversional canal systems (Pipe systems)

% k = loop for outlet locs
% l = loop for inlet elevation levels
% m = loop for inlet per elevation level
% ndl = loop for Q decreaser
% j = loop for dam heights

% Relevant parameters for both plant types
cellsz_m=500;
minQ_small = 0.1;           %Minimum design discharge (in m3/s) required in inlet for smallDP
minQ_large = 1.0;           %Minimum design discharge (in m3/s) required for inlet in largeDP and outlet in large RP
minQ=0.1;                   %Minimum annual avg flow to be valid inlet (in m3/s)
dowin=50;                   %Winsize for DP and P calculations.
do=do_main_km*1e3/cellsz_m; %Distance between the outlets (spacing), in terms of number of cells. 1 cell is ~500m. Standard:50 (25km); for development:2000(1000km)
inlet_main_sradius_km=(do_main_km-minDis_main_km);
inlet_trib_sradius_km=(do_trib_km-minDis_trib_km);
inlet_main_sradius=inlet_main_sradius_km*1e3/cellsz_m;   %Search radius tributary P systems in # of cells
inlet_trib_sradius=inlet_trib_sradius_km*1e3/cellsz_m;   %Search radius main P systems in # of cells
cost_lim=0.5;               %Cost limit $/kWh
outletdeselect=2;           %Number of deselected sites upstream of outlet

%buffers
buffer_cult_km = 1;         %Buffer for cultural heritages in km
buffer_nat_km = 2;          %Buffer for natural heritages in km
buffer_thrust_km = 1;       %Buffer for thrust in km
buffer_fault_km = buffer_thrust_km;        %Buffer for fault in km
buffer_glof_km = 0.5;       %Buffer for glof paths in km

% Costfactors
ResettleFactor=3;           %GDP scale factor for per capita resettlement cost
CommBenefitRate_RP=10/100;     %Community benefit sharing or compensation added as % of resettlement cost
CommBenefitRate_DP=10/100;
LandAcqRate_RP=15/100;      %Added price for land acquisition by govt
LandAcqRate_DP=15/100;      %Added price for land acquisition by private sector
geohazard_level_rates =[0	0  % Added cost for each hazard level 1=Very Low, 5=Very High
    1	2
    2	4
    3	6
    4	8
    5	10];

% RD parameters
%RD_minH=4;                  %Minimum dam height or head for RD project in m
RDwinsize_hires_dem=30;     %Search radius for dam width in number of cells in finer DEM ~ 3km = 3000m/100m=30 cells
HRcellsz_m=100;             %Resolution of high res DEM in m
%storage_lim=5/100;          %Maximum limit on reservoir size as percentage of annual average inflow volume

% P parameters
nd=4;                       %Number of Q-decreaser loops
n_ielevations=6;            %Number of inlet elevation levels per outlet loc
% Plarge_minH=4;              %Minimum head in m for large P (default =4m)
% Psmall_minH=20;             %Minimum head in m for small P (default =20m)

% Other params
rivermouth_inland=400;      %Number of cells inland of basin outlet deselected (200km) (rivermouth_constraint)
depth_cutoff=4;             %Water depth cutoff (m)
bighydro_cutoff=50;         %Max Q based on small-scale Hydropower Veileder 10MW document p20
MiniDamMax=30;              %Max Dam Height based on Q cutoff. This is consistent with big Hydro Veileder p62

%Disabled many of David's constraints that are not relevant
% Constraints. 0=constraint off / 1= constraint on
runDP=0;                    %Run dam-pipe-resevoir systems 0=no/1=yes
runP=1;                     %Run pipe systems 0=no/1=yes
runRD=1;                    %Run river dam system 0=no/1=yes
runMini=0;                  %Run mini dam systems 0=no/1=yes

% Tech constraints - yes(1) or no(0)
%enablereslimit = 1;              %Restrict reservoir storage volume to be a % of annual average inflow volume
varyRPDP_Qdesign = 1;                %Use different Q design for RD and DP project types: (0=) use first val(1=) use all vals in sel_design_flows
waterbodies_constr=1;           %No dams in existing water bodies
rivermouth_constr=0;        %No dams in rivermouth (see outletdeselect-setting)
mainstream_constr=0;        %Deselect dam locations on basin mainstream
deselect_mainstream=0;      %Deselect DP and P locations that are on the mainstream of RD power station.
navi_constr=0;              %No dams in navigable rivers (see depth_cutoff-setting)
BeforeFirstDam_constraint=0;%No dams on mainstream locations before the first GrandDam
MiniHydro_deselect=0;       %Deselect dam locations on small streams (MiniHydro locations) (see bighydro_cutoff-setting)
MiniHydro_select=0;         %Select dam locations on small streams (MiniHydro locations) (see bighydro_cutoff-setting)
MiniHydro_special=0;        %Special costmodel activation for MiniHydro locations (see bighydro_cutoff-setting and MiniDamMax-setting)
mangrove_constr=0;          %No dams in mangroves

% Econ constraints - yes(1) or no(0)
cost_constr=0;          %Enabling max cost constraint (see cost_lim-setting)
discost_set=1;          %Taking into account distance to load cost. 0=no/1=yes
landval_set=1;          %Taking into account Land value cost. 0=no/1=yes

ASYMRD=0;                   %Asymmetric dams RD-systems 0=no/1=yes
ASYMDP=0;                   %Asymmetric dams DP-systems 0=no/1=yes
betaC=0;                    %Deselection cost weight
betaP=0;                    %Deselection power weight
betaCP=1;                   %Deselection cost per power weight
betaA=0;                    %Deselection flow acc (as indication of river lenght)
betaL=0;                    %Deselection based on Lake size. Cost per power is always on. This setting should be on only in eco scenario

%% Setup output runnames and filenames and start log
% runname =<prefix>_<scenariotype>_<policytype>_<constraintype>_<hazard represntation type>
if runsustain
    if ~all(geohazard_select)
        runname_prefix=strjoin([runname_prefix(:)', geohazard_snames(find(geohazard_select))],'_'); %update to have disaster in name
    end
    runname= strjoin([runname_prefix,scenario,policynames{policytype+1},"Sust",geohazard_scennames{geohazard_scentype+1}],'_');
else
    runname=strjoin([runname_prefix,scenario, policynames{policytype+1}, "Tech_Fin"],'_');
end

if savedata
    matpath = fullfile(rootpath, 'output', scenario, continent_in, runname);
    if ~isfolder(matpath); mkdir(matpath); end
    logf = fullfile(matpath, sprintf('Hydrus_%s_%s.log',string(datetime('now','Format','dd-MMM-yyyy')),runname));
    diary(logf)
    % Setup output file locations
    matfile = fullfile(matpath, sprintf('Basin%d_output_do.mat', nbasin));
    matfileCOEPOT = fullfile(matpath, sprintf('COEPOT_b%d_do.mat', nbasin));
    xlsfile=strrep(matfileCOEPOT,'mat','xlsx');
end

return %temp for testing parfor loops

%% Read data file for current basin
basinfname = fullfile(rootpath, 'data', continent_in, 'Basin_UIB', sprintf('PantpeBasin_%d.mat', nbasin));
load(basinfname);
coordsys=Rw.CoordinateSystemType; % Get coordinate system from data
[nrw,ncw] = size(acc);
disp("Input data loaded!!!")

%% Display start messages
mystartstatement='HYDRUS02 START RUN for %s: \n %s\n';
fprintf(mystartstatement,runname,datetime('now'));
%
disp('Start')
fprintf('\nTECH-ECON SETTINGS\n');
fprintf('Continent: %s\n', continent_in);
fprintf('Basin: %d \n', nbasin);
fprintf('Coordinate system: %s \n',coordsys);
fprintf('Scenario: %s\n', scenario);
fprintf('Policy scenario type: %s\n', policynames{policytype+1});
fprintf('Vary plant types, do and sradius based on # of river levels: %d\n', varyRiverLevel_HPtype);
fprintf('Vary Q design between R and P: %d\n', varyRPDP_Qdesign);
fprintf('Vary Q design between small and large P: %d\n', enablesmallPdesignQH);
fprintf('Tributaries Spacing distance in km: %d\n', do_trib_km);
fprintf('Mainstream Spacing distance in km: %d\n', do_main_km);
fprintf('Tributaries Sradius in km: %0.1f\n', inlet_trib_sradius_km);
fprintf('Mainstream Sradius in km: %d\n', inlet_main_sradius_km);
fprintf('Existing dams constraint: %d\n', ExistingDams_constraint);
fprintf('\nSUSTAIN SETTINGS\n');
fprintf('Discharge corrected for water consumption: %d\n', waterconsumption);
fprintf('Slack flow constraint: %d%%\n', slackflow_constraint*slackflow*100);
fprintf('Protected areas constraint: %d\n', protarea_constr);
fprintf('Geohazard risk scenario type: %s\n', geohazard_scennames{geohazard_scentype+1});
fprintf('Geohazards selected: \n')
disp( geohazard_names(find(geohazard_select)))
% fprintf('Navigational constraint: %d\n', navi_constr);
% fprintf('No dams before first dam: %d\n', BeforeFirstDam_constraint);
%fprintf('Deselection based on lake size: %d\n', betaL);

% Reading ancillary scripts
disp('Reading ancillary scripts');
addpathname1 = fullfile(rootpath,'Hydrus','functions');
addpathname2 = fullfile(rootpath,'Hydrus','scripts');
addpath(addpathname1, addpathname2)

% add routing script
defdirs
% some kind of check point
returnID=0;


%% Update LandVal and Protected Area datasets
LandValue = LandAgriVal + LandAcqVal + LandTreeVal;
ProtArea = createBuffer(Heritages_natural, buffer_nat_km*1e3/cellsz_m)+ ...
    createBuffer(Heritages_cultural, buffer_cult_km*1e3/cellsz_m);

%% Geo-hazard: Setup geo-hazard scenario settings
if runsustain
    if geohazard_scentype==0  % No geohazard
        geohazard_cost      =0;
        geohazard_constr    =0;
    elseif geohazard_scentype==1 % Status Quo
        geohazard_cost      =1; %Add costs for hazard levels
        geohazard_constr    =0; %Skip areas with unaccetable hazard levels: (1 Risk Averse scen=) skip high max hazard areas; (2 multi hazard=) skip high mult hazard areas
    elseif geohazard_scentype==2 % Risk averse  or binary
        geohazard_cost      =1; %Add costs for hazard levels
        geohazard_constr    =1; %Skip areas with unaccetable hazard levels: (1 Risk Averse scen=)skip high max hazard areas; (2 multi hazard=) skip high mult hazard areas
    elseif geohazard_scentype ==3 %Multi hazard
        geohazard_cost      =1; %Add costs for hazard levels
        geohazard_constr    =2; %Skip areas with unaccetable hazard levels: (1 Risk Averse scen=)skip high max hazard areas; (2 multi hazard=) skip high mult hazard areas
    end
end 
%% Geo-hazard: Get hazard cost rate for different geohazard scenes
% Compile hazard risk data
geohazard_compile(:,:,1)=geohazard_select(1)*seismic_level;           % seismic_5level
geohazard_compile(:,:,2)=geohazard_select(2)*createBuffer(seismicthrust>0, buffer_thrust_km*1e3/cellsz_m)*nhazardlevel;  % seismicthrust_level considered highest risk level
geohazard_compile(:,:,3)=geohazard_select(3)*landslide_level;          % landslide_level Reclassified
geohazard_compile(:,:,4)=geohazard_select(4)*createBuffer(glofpath,buffer_glof_km*1e3/cellsz_m)*nhazardlevel;  % glof_levelGLOF  path + buffer considered highest risk level

if geohazard_cost && geohazard_constr<=1 % for geohazard status quo and risk averse scenarios
    % Based on max hazard class
    % Composite considering max of all hazards
    geohazard_max_levels=maskBasin(max(geohazard_compile,[],3), ~outside);
    
    % Get rates for hazard classes
    geohazard_hazard_rates=changem(geohazard_max_levels,geohazard_level_rates(:,2),geohazard_level_rates(:,1));
    
elseif geohazard_cost && geohazard_constr ==2 %% for multi-hazard scenario
    % Composite considering max of all hazards rescaled from 0-5
    geohazard_sum_naned=maskBasin(sum(geohazard_compile,3), ~outside);
    
    geohazard_multi_score=rescale(geohazard_sum_naned,0,5,'InputMin',0,'InputMax',20);
    geohazard_multi_levels=reclassify(geohazard_multi_score, 0:4,0:5); %,hazard_cnames);
    %Get rates for hazard classes
    geohazard_hazard_rates=changem(geohazard_multi_levels,geohazard_level_rates(:,2),geohazard_level_rates(:,1));
end

%% Geo-hazard: Get geohazard areas to skip
if geohazard_constr==1 % Risk averse case
    geohazard_skip=geohazard_max_levels==5; % skip only very high level
elseif geohazard_constr==2 % Multi-hazard case
    geohazard_skip=geohazard_multi_levels>=4; % skip high and very high level
end

%% Setup selected Qdesign Q-maps with or without water consumption
if waterconsumption
    Qsuffix='_wc';
    Qyravg=Qwc;
else
    Qyravg=Q;
    Qsuffix='';
end

% Select RD design flow
lvl=1;
Q_RD_design_sel=eval(strcat(sel_design_flows{lvl}, 'design', Qsuffix));
Q_RD_design_sel_LF=eval(strcat(sel_design_flows{lvl}, 'design_LF', Qsuffix));

% Select largeDP design flow
if varyRPDP_Qdesign
    % Get Q for Ptype projects, give whole matrix as necessary for inlet search
    lvl=2;
    Q_P_design = eval(strcat(sel_design_flows{lvl}, 'design', Qsuffix));
    Q_P_design_LF = eval(strcat(sel_design_flows{lvl}, 'design_LF', Qsuffix));
    Q_P_design_mean = eval(strcat(sel_design_flows{lvl}, 'design_mean', Qsuffix));
else % use same Q as RD
    Q_P_design = Q_RD_design_sel;
    Q_P_design_LF = Q_RD_design_sel_LF;
    Q_P_design_mean = eval(strcat(sel_design_flows{lvl}, 'design_mean', Qsuffix));
end

% Select smallDP design flow
if enablesmallPdesignQH
    lvl=3;
    Q_P_design_small = eval(strcat(sel_design_flows{lvl}, 'design', Qsuffix));
    Q_P_design_LF_small = eval(strcat(sel_design_flows{lvl}, 'design_LF', Qsuffix));
    Q_P_design_mean_small = eval(strcat(sel_design_flows{lvl}, 'design_mean', Qsuffix));
end
fprintf("Qdesign loaded as: \n\t\tRP: %s \n\t\tlargeDP: %s \n\t\tSmallDP: %s\n",sel_design_flows{:})

%% Generate outlets and their IDs
disp('Assigning powerstation locations');
%%%%find_outlets_wins.m START

% Skip if max flowdist is < do (distance between outlets)
if max(flowdist(:)) < do
    disp('Flowdist is less than distance between outlets');
    save(matfile,'-v7.3','do');
    returnID=1;
    return
end

if varyRiverLevel_HPtype==0
    [ro,co] = find_outletsInChannel(channel_main_trib>0,flowdist,do_main_km*1e3/cellsz_m);
    out_type = 2*ones(size(ro)); %all=2
    out_type_RD=2; %RP and DP made in both
elseif varyRiverLevel_HPtype==1
    [r1,c1] = find_outletsInChannel(channel_main_trib==1,flowdist,do_trib_km*1e3/cellsz_m); %tributaries
    [r2,c2] = find_outletsInChannel(channel_main_trib==2,flowdist,do_main_km*1e3/cellsz_m); %mainstreams
    ro=[r1;r2];
    co=[c1;c2];
    out_type = [ones(size(r1)); 2*ones(size(r2))]; %tributatries=1, mainstream=2;
    out_type_RD=2; %RP only in mainstreams
elseif varyRiverLevel_HPtype==2
    [r1,c1] = find_outletsInChannel(channel_ord==1,flowdist,do_trib_km*1e3/cellsz_m); %tertiary
    [r2,c2] = find_outletsInChannel(channel_ord==2,flowdist,do_main_km*1e3/cellsz_m); %secondary
    [r3,c3] = find_outletsInChannel(channel_ord==3,flowdist,do_main_km*1e3/cellsz_m); %primary
    ro=[r1;r2;r3];
    co=[c1;c2;c3];
    out_type = [ones(size(r1)); 2*ones(size(r2)); 3*ones(size(r3))]; %tertiary=1, secondary=2, primary=3
    out_type_RD=3; %RP only in primary
end
outlets = sub2ind(size(acc),ro,co);

% archive initial list for plotting in the end
ro_arc=ro;
co_arc=co;

%% Make variable placeholders
noutlets=numel(ro);
Accoutlets=zeros(1,noutlets);
PopDisplacedOpt=zeros(1,noutlets);
hazard_rate_RD=zeros(1,noutlets);
COETotRD=zeros(1,noutlets);
RDPnet=zeros(1,noutlets);
RDtheory_GWh=zeros(1,noutlets);
Q_RD_design=zeros(1,noutlets);
Q_RD_design_LF=zeros(1,noutlets);
RDRegion_id=zeros(1,noutlets);
RDCountry_id=zeros(1,noutlets);
latOut=zeros(1,noutlets);
lonOut=zeros(1,noutlets);
xOut=zeros(1,noutlets);
yOut=zeros(1,noutlets);
OptDH=zeros(1,noutlets);
OptDL=zeros(1,noutlets);
OptPop=zeros(1,noutlets);
OptLV=zeros(1,noutlets);
RDDepth=zeros(1,noutlets);
OptSpecCap=zeros(1,noutlets);
OptInv=zeros(1,noutlets);
RDP=zeros(1,noutlets);
RDlakeSurface=cell(1,noutlets);
RDVolumeLake=cell(1,noutlets);
RDVolumeLake15s=zeros(1,noutlets);
RDSurfaceLake15s=zeros(1,noutlets);
Zoutlets=zeros(1,noutlets);
DisTransOutlet=zeros(1,noutlets);
DisRoadOutlet=zeros(1,noutlets);
DisSettlementOutlet=zeros(1,noutlets);
ZiWinfr=zeros(1,noutlets);

%% Get coordinates (Latlon or XY) for outlets
for  k = 1:numel(ro)
    if strcmp(coordsys, 'geographic')
        [latOut(k), lonOut(k)] = setltln(acc, Rw, ro(k), co(k)); %Coordinates outlets
    elseif strcmp(coordsys, 'planar')
        % convert r,c to world xy and then lat lon
        [xOut(k), yOut(k)] =  intrinsicToWorld(Rw,co(k), ro(k)) ;
        [latOut(k), lonOut(k)] = projinv(proj,xOut(k), yOut(k)); %Coordinates of outlets
    end
    Depth_Rivermouth(k) = D(ro(k),co(k)); % Depth of locations from outlet
    flowdist_rivermouth(k) = flowdist(ro(k),co(k)); % Distance from outlet
    %Qoutlets(k) = Qyravg(ro(k),co(k));
    outlets(k) = sub2ind(size(Z),ro(k),co(k));
    Accoutlets(k) = acc(ro(k),co(k));
end

%% Apply binary constraints - Tech/Fin
% Deselect outlets with Qyravg<= minQ - this not necessary anymore actually
%  coz channel already only has Qannual_m3s>=0.1
for k= 1:numel(ro)
    if isnan(ro(k))==1; continue; end;
    if Qyravg(ro(k),co(k))< minQ
        ro(k)=NaN;
        co(k)=NaN;
    end
end

% Deselect outlets that fall on existing water bodies
if waterbodies_constr ==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        if AllWaterBodies(ro(k),co(k))== 1
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
end

% Rivermouth constraint - Deselect outlet and - No dam in river mouth
if rivermouth_constr==1
    if numel(ro)<outletdeselect
        disp('To few ro and co');
        save(matfile,'-v7.3','do');
        returnID=1;
        return
    end
    [~,sortQidx]=sort(Accoutlets,'descend');
    
    %First two for sure
    ro(sortQidx(1:outletdeselect))=NaN;
    co(sortQidx(1:outletdeselect))=NaN;
    %     ro(sortQidx(50:(end)))=NaN;
    %     co(sortQidx(50:(end)))=NaN;
    
    %This method causes additional deselection after lakes in the river
    %system, especially in Congo and Amazon
    %Second deselect dams on mainstream xkm inland and on rivers with a certain depth
    for k=1:numel(ro)
        if flowdist_rivermouth(k) < rivermouth_inland && Depth_Rivermouth(k) > depth_cutoff
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
    
    %Other 2 maybe depending on river depth
    %     if strcmp('SAM',continent_in)
    %         if nbasin==1
    %             if numel(sortQidx) < 8; RDi = numel(sortQidx); else RDi=33; end
    %             for i=3:RDi
    %                 if Depth_Rivermouth(sortQidx(i)) > depth_cutoff
    %                     ro(sortQidx(i))=NaN;
    %                     co(sortQidx(i))=NaN;
    %                 end
    %             end
    %         end
    %     end
end
% sum(isnan(ro))

% figure(1); hold on
% plot(co,ro,'b.','markersize',15); hold off
% figure(2); subplot(1,3,2); hold on
% plot(co,ro,'b.','markersize',15); hold off
% fd=flowdist;fd(Q>100)=0;cmap=jet(256);cmap(1,:)=[1 1 1];
% figure(2);ax3=subplot(1,3,3); imagesc(fd);colormap(cmap);axis image;

% Navigation constraint for navigability (4m)
if navi_constr==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        if D(ro(k),co(k))> depth_cutoff
            ro(k)=NaN;
            co(k)=NaN;
        else
        end
    end
else
end
% [rd, cd] = setpostn(Q, Rw, -3.128930, -60.027996);
% figure(2); hold on
% plot(cd,rd,'g.','markersize',25);
% plot(co,ro,'b.','markersize',15); hold off
% sum(isnan(ro))

% Mangrove constraint
if mangrove_constr==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end
        if MC(ro(k),co(k))==1
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
end

% Mainstream-constraint - No dams on basin mainstream
if mainstream_constr==1
    [~,MainIdx] = max(Qyravg(:));
    mainstrm = find_mainstream(Qyravg,MainIdx,adir);
    outlets(ismember(outlets,mainstrm))= NaN;
    ro(find(isnan(outlets)))=NaN;
    co(find(isnan(outlets)))=NaN;
end

% MiniHydro deselection - Disable small plants or outlets w Q< bighydro_cutoff
if MiniHydro_deselect==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        if Qyravg(ro(k),co(k))< bighydro_cutoff
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
end

% MiniHydro Selection - Disable big plants or outlets w Q> bighydro_cutoff
if MiniHydro_select==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        if Qyravg(ro(k),co(k))> bighydro_cutoff
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
end

%% Existing dams de-selection
if ExistingDams_constraint==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        % Disable outlets where existing dams lie
        if NoDamsLand(ro(k),co(k))
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
else
    % Disable GrandDatabase location for full scenario or
    % ExistingDams_constraint==0
    %nullify existing dams database because overlap search code calls these later
    r_damst=-99; c_damst=-99; GrandIdx=-99; Grandlat=-99; Grandlon=-99;
    NoDamsLand=0*NoDamsLand;
end

% No dams before first GrandDam on mainstream
if BeforeFirstDam_constraint==1
    if r_damst(1)~=-99;
        NoBeforeFirstMainstrDam = FirstMainstrmDam(Qyravg,adir,acc,r_damst,c_damst,GrandIdx);
        outlets(ismember(outlets,NoBeforeFirstMainstrDam))= NaN;
        ro(find(isnan(outlets)))=NaN;
        co(find(isnan(outlets)))=NaN;
    end
end

% if all NaN break loop
if sum(isnan(ro))==numel(ro)
    disp('First all NaN break');
    save(matfile,'-v7.3','do');
    returnID=1;
    return
end

%% Apply binary constraints - Sustainable
% Deselect locations in protected areas (natural + cultural)
if protarea_constr==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end
        if ProtArea(ro(k),co(k)) % for all non 0 vals
            ro(k)=NaN;
            co(k)=NaN;
        else
        end
    end
else
end

% Geohazard-risk constraint - No plants in selected geohazard risk areas
if geohazard_constr %can be 1 or 2 for risk averse and multi-hazard scenario cases
    for k= 1:numel(ro)
        if isnan(ro(k)); continue; end
        if geohazard_skip(ro(k),co(k))  % for all non 0 vals
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
end

% Geohazard-cost constraint - cost increase for medium hazard areas
if geohazard_cost==1
    for k=1:numel(ro)
        if isnan(ro(k)); hazard_rate_RD(k)=NaN; continue; end
        if ~isnan(geohazard_hazard_rates(ro(k),co(k)))
            hazard_rate_RD(k) = geohazard_hazard_rates(ro(k),co(k)); % Added cost associated w hazard risk
        end
    end
else
    hazard_rate_RD(:) = 0;
end

%% Extract relevant data values at valid outlets
for k=1:numel(ro)
    if isnan(ro(k)); continue; end
    
    Zoutlets(k) = Z(ro(k),co(k));
    RDRegion_id(k) = Regions(ro(k),co(k));
    RDCountry_id(k) = Countries(ro(k),co(k));
    RDDepth(k) = D(ro(k),co(k));
    DisTransOutlet(k) = DisTransmission(ro(k),co(k)); % Distance to powerline (km)
    DisRoadOutlet(k) = DisRoad(ro(k),co(k)); % Distance to road (km)
    DisSettlementOutlet(k) = DisSettlement(ro(k),co(k)); % Distance to nearest settlement (km)
    
    % Get Q for RDtype project analysis at each outlet
    Q_RD_design(k) = Q_RD_design_sel(ro(k),co(k));
    Q_RD_design_LF(k) = Q_RD_design_sel_LF(ro(k),co(k));
end

%% Slackflow for ecological reasons - not applied on smallP
if slackflow_constraint==1
    Q_RD_design = Q_RD_design *( 1-slackflow);
    Q_P_design = Q_P_design *( 1-slackflow);
    Q_P_design_mean = Q_P_design_mean *( 1-slackflow);
end

%% Report deselection of outlets
n_outlets=numel(outlets);
n_outlets_1stdeselect=sum(isnan(ro));
fprintf('Number of outlets: %d.\n', n_outlets);
fprintf('Number of deselected: %d.\n', n_outlets_1stdeselect);


% figure(labindex +1);clf;imagescnan(flowdist);colormap(flipud(gray));axis image; hold on
% plot(co_arc,ro_arc,'ko','markersize',10,'DisplayName','Outlets: All');
% hold on
% plot(co,ro,'g.','markersize',15,'DisplayName','Outlets: Selected')
% legend('-DynamicLegend','Location','best')
% title(sprintf('%d outlets deselected out of %d \n',n_outlets_1stdeselect,n_outlets))
% hold off

% %% Control for NaNs in Dis map (slightly inconsistent with other maps because different sea map was used)
% for k=1:numel(DisTransOutlet)
%     if isnan(DisTransOutlet(k))==1; ro(k)=NaN; co(k)=NaN; end;
% end

% %% Control for NaNs in Dis map (slightly inconsistent with other maps because different sea map was used)
% for k=1:numel(DisTransOutlet)
%     if isnan(DisTransOutlet(k))==1; ro(k)=NaN; co(k)=NaN; end;
% end

%% Create windows for Qdecreaser
for k = 1:numel(outlets)
    if isnan(ro(k)); continue; end
    
    inlet_win = dowin;
    ZiWinfr(k) = ro(k)-inlet_win;
    ZiWinlr(k) = ro(k)+inlet_win;
    ZiWinfc(k) = co(k)-inlet_win;
    ZiWinlc(k) = co(k)+inlet_win;
    
    if ZiWinfr(k)<1 || ZiWinlr(k)>nrw || ZiWinfc(k)<1 || ZiWinlc(k)>ncw
        ro(k)=NaN;
        co(k)=NaN;
    else
    end
    
    if isnan(ro(k)); continue; end
    
    fdir_inlet_win{k}      = fdir(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    adir_inlet_win{k}      = adir(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    Q_inlet_win{k}         = Qyravg(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    acc_inlet_win{k}       = acc(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    Z_inlet_win{k}         = Z(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    %Pop_inlet_win{k}       = Pop(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    protarea_inlet_win{k}  = ProtArea(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    WDPA_PL20_inlet_win{k} = 0*protarea_inlet_win{k};
    %     LandValue_inlet_win{k} = LandValue(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    %     Dis_inlet_win{k}       = DisTransmission(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    flowdist_inlet_win{k}  = flowdist(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    
end

WinProb=sum(isnan(ro))- n_outlets_1stdeselect;
fprintf('Number of window problems: %d.\n', WinProb);

index_inlet_win = sub2ind([(2*inlet_win+1) (2*inlet_win+1)],inlet_win+1,inlet_win+1);

% initialize variable
n_DPinlets=zeros(1,k); %number of inlets tried in each outlet over all elevations

%% If no valid outlets break loop
if sum(isnan(ro))==numel(ro)
    disp('all NaN brak loop');
    save(matfile,'-v7.3','do');
    save(matfileCOEPOT,'-v7.3','do');
    returnID=1;
    return
end

%% NOT RUN RIGHT NOW -- DP: RESERVOIR-SYSTEMS, WITH DAM AND PIPE, LOADFACTOR 100% ---  as runDP=0
% if runDP==1
%     for ndl = 1:nd
%
%         fprintf('Q decreaser progress: %d of %d.\n', ndl, nd);
%
%         %% Find inlets
%         disp('Find reservoir locations');
%
%         find_inlets
%
%         %% Determine lat lon coordinates of outlets and inlets
%         disp('Calculate pipelenght');
%
%         pipelength
%
%         %% Determine lat lon coordinates of outlets and inlets
%         disp('Find lat lon coordinates');
%
%         latlon
%
%         %% Calculate dam height
%         disp('Calculate dam height and width');
%
%         dam_height_width
%
%         %% Calculate Lake
%         disp('Determine lake');
%
%         DPlake
%
%         %% Running cost model
%         disp('Running cost model');
%
%         damcostmodel
%
%         %% Selecting best inlet location
%         disp('Selecting best reservoir-inlet location');
%
%         inletselecter
%
%         %% Q-decreaser
%         disp('Q decreaser');
%
%         Qdecreaser
%
%         %% Store info per Q decreaser loop
%         disp('Store info per Q decreaser loop')
%
%         COEend{ndl}=aCOEmin; % $/kWh Dam-Pipe COE
%         DPPnetend{ndl}=aPnetmin; % GWh Dam-Pipe energy potential
%         latminend{ndl} = lat_in_min;
%         lonminend{ndl}= lon_in_min;
%         rinMinend{ndl} = rinMin;
%         cinMinend{ndl}= cinMin;
%         aInletminEnd{ndl}=aInletMin;
%         ainlet_windowMinEnd{ndl}=ainlet_windowMin;
%         ahdammin{ndl}=hdammin;
%         aldammin{ndl}=ldammin;
%         DPDepth{ndl}=DPDepth2;
%         CostElementsMinend{nd}=CostElementsMin;
%
%         aDPPnetEnd3  = horzcat(DPPnetend{:})'; % GWh  Dam-Pipe systems Pnet
%         fprintf('Total DP potential: %0.0f GWh\n',sum(aDPPnetEnd3(~isnan(aDPPnetEnd3))));
%
%     end
% end

%% P: RUN-OF-RIVER PIPE-SYSTEMS, NO DAM, LOADFACTOR <100%
if runP==1
    for ndl = 1:nd
        
        fprintf('\n\n\n>>>>>>>Pipe Systems Q decreaser: %d of %d.\n', ndl, nd);
        
        %%%% PipeSystems.m START
        %% Pipe-systems procedure
        
        %% Find pipe inlets
        disp('Find Pinlets');
        
        %%%%find_Pinlets.m START
        %% Find Pinlet locations
        for k = 1:numel(outlets)
            if isnan(ro(k)); Pinlet_win{k}=[]; continue; end
            
            %Data windows are created in find_outlet_wins
            
            [nr_ZiW,nc_ZiW] = size(Z_inlet_win{k});
            
            Z_inlet_win_temp=0;
            rr1=inlet_win+1-inlet_trib_sradius;
            rr2=inlet_win+1+inlet_trib_sradius;
            cc1=inlet_win+1-inlet_trib_sradius;
            cc2=inlet_win+1+inlet_trib_sradius;
            Z_inlet_win_temp = Z_inlet_win{k}(rr1:rr2,cc1:cc2);
            
            Zmin = Z_inlet_win{k}(inlet_win+1,inlet_win+1);
            Zmax = max(Z_inlet_win_temp(:));
            Zupdown = Zmax-Zmin;
            % Get 6 elevation levels within the min/max of the data window
            ZPin = linspace(Zmin+(0.1*Zupdown),Zmax-(0.1*Zupdown),n_ielevations);
            
            Poutlet_win=  sub2ind(size(Z_inlet_win{k}),inlet_win+1,inlet_win+1); %index of center dam outlet location
            %Zout=  sub2ind(size(Z_inlet_win{k}),inlet_win+1,inlet_win+1); %index of center dam outlet location
            
            %figure(2);clf;imagesc(Z_inlet_win{k});axis image; colormap(flipud(gray));hold on
            %plot(inlet_win+1,inlet_win+1,'.r','markersize',20); hold off
            
            for l = 1:numel(ZPin)
                % vary search radius based on out_type if split
                if varyRiverLevel_HPtype && out_type(k)==1   %in tributaries when split
                    sel_inlet_sradius=inlet_trib_sradius;
                else  %if main or dont split then use inlet_main_sradius
                    sel_inlet_sradius=inlet_main_sradius;
                end
                Pinletsw = inlet_finderSD(Z_inlet_win{k},Q_inlet_win{k},acc_inlet_win{k},...
                    fdir_inlet_win{k},protarea_inlet_win{k},WDPA_PL20_inlet_win{k},protarea_constr,...
                    nr_ZiW,nc_ZiW,Poutlet_win,ZPin(l),deselect_mainstream,index_inlet_win,adir_inlet_win{k},...
                    flowdist_inlet_win{k},sel_inlet_sradius);
                
                Pinlet_win{k}{l} = find(Pinletsw);
                %convert inlet maps to arrays of rows/columns
                rPin_w{k}{l} = rem(Pinlet_win{k}{l}-1, nr_ZiW)+1;
                cPin_w{k}{l} = fix((Pinlet_win{k}{l}-1)/nr_ZiW)+1;
                % Converting window r,c to full catchment r,c
                rPin{k}{l} = ro(k)+(rPin_w{k}{l}-inlet_win-1);
                cPin{k}{l} = co(k)+(cPin_w{k}{l}-inlet_win-1);
                Pinlet{k}{l} = sub2ind(size(Z),rPin{k}{l},cPin{k}{l});
                ZPinlet{k}{l} = Z(Pinlet{k}{l});
                QPinlet{k}{l} = Qyravg(Pinlet{k}{l});
                accPinlet{k}{l} = acc(Pinlet{k}{l});
                PRegion_id{k}{l} = Regions(rPin{k}{l},cPin{k}{l});
                PCountry_id{k}{l} = Countries(rPin{k}{l},cPin{k}{l});
                
                %% Extract Q at inlet
                %Set Qtrib to Qsmall if trib outlet and enablesmallP is on,
                % for main outlet or when enablesmall is off use Q30
                if enablesmallPdesignQH==1 && out_type(k)==1
                    QDesignPinlet{k}{l} = Q_P_design_small(Pinlet{k}{l});
                    QDesignLFPinlet{k}{l} = Q_P_design_LF_small(Pinlet{k}{l});
                    QDesignMeanPinlet{k}{l} = Q_P_design_mean_small(Pinlet{k}{l});
                else
                    QDesignPinlet{k}{l} = Q_P_design(Pinlet{k}{l});
                    QDesignLFPinlet{k}{l} = Q_P_design_LF(Pinlet{k}{l});
                    QDesignMeanPinlet{k}{l} = Q_P_design_mean(Pinlet{k}{l});
                end
                %
                %                 if slackflow_constraint==1
                %                     QDesignPinlet{k}{l} = QDesignPinlet{k}{l} *( 1-slackflow);
                %                 end
                %
                %Check
                %figure(3);clf;imagesc(log(Q_inlet_win{3}));axis image; colormap(flipud(gray));hold on
                %plot(inlet_win+1, inlet_win+1,'b.','markersize',24);
                %plot(cPin_w{3}{1}, rPin_w{3}{1},'r.','markersize',24);
                %hold off
                
                %figure(3);clf;imagesc(log(Q));axis image; colormap(flipud(gray));hold on
                %plot(cPin{3}{1}, rPin{3}{1},'r.','markersize',24);
                %hold off
                
            end
        end
        
        %% Variable creator
        for k = 1:numel(outlets)
            
            PopDP{k}=[];
            PopCostDP{k}=[];
            LandValueDP{k}=[];
            COEP{k}=[];
            PL{k}=[];
            latP{k}=[];
            lonP{k}=[];
            PP_energy_GWh{k}=[];
            PP_Theory_GWh{k}=[];
            CostElementsP{k}=[];
            nPipeP{k}=[];
            PPipeDia{k}=[];
            OptInvP{k}=[];
            PP_W{k}=[];
            bMinElevP(k)=0;
            aCOEPmin(k)=NaN;
            aPPnetmin(k)=NaN;
            aPPtheorymin(k)=NaN;
            nPipeminP(k)=NaN;
            CostElementsPMin{k}=[];
            dfQPmin(k)=NaN;
            accPmin(k)=NaN;
            dfPLmin(k)=NaN;
            HeadPmin(k)=NaN;
            HeadraceLPmin(k)=NaN;
            aPInletMin(k)=NaN;
            aPinlet_windowMin(k)=NaN;
            lat_Pin_min(k)=NaN;
            lon_Pin_min(k)=NaN;
            rPinMin(k)=NaN;
            cPinMin(k)=NaN;
            QDesignPinletMin(k) =NaN;
            QDesignLFPinletMin(k) =NaN;
            QDesignMeanPinletMin(k) = NaN;
            ZPinletMin(k) =NaN;
            nPipePMin(k) =NaN;
            OptInvPMin(k) =NaN;
            PPMin(k) =NaN;
            DeselectedSites_unSortRD{k}=[];
            DeselectedGrandSites_unSortRD{k}=[];
            MiniFlag(k) = 0;
            OptInv(k) =NaN;
            OptSpecCap(k) =NaN;
            RDP(k) = NaN;
            RDVolumeLake15s(k) = NaN;
            RDSurfaceLake15s(k) = NaN;
            OptPop(k) = NaN;
            OptLV(k) = NaN;
            OptSpecCapP{k} = [];
            OpthfP{k} = [];
            OptDP{k} = [];
            OptSpecCapPMin(k) = NaN;
            OpthfPMin(k) = NaN;
            PSurfaceLake15s(k) = 1e-3;
            OptDPMin(k) = NaN;
            
            if numel(Pinlet_win{k})==0; continue; end
            
            for l = 1:n_ielevations
                COEP{k}{l}=[];
                PL{k}{l}=[];
                latP{k}{l}=[];
                lonP{k}{l}=[];
                PP_energy_GWh{k}{l}=[];
                PP_Theory_GWh{k}{l}=[];
                PPipeDia{k}{l}=[];
                OptInvP{k}{l}=[];
                PP_W{k}{l}=[];
                COEPminElev{k}(l)=0;
                bMinP{k}(l)=0;
                OptSpecCapP{k}{l}=[];
                OpthfP{k}{l}=[];
                OptDP{k}{l}=[];
                
                for m = 1:numel(Pinlet_win{k}{l})
                    COEP{k}{l}(m)=NaN;
                    PL{k}{l}(m)=NaN;
                    latP{k}{l}(m)=NaN;
                    lonP{k}{l}(m)=NaN;
                    PP_energy_GWh{k}{l}(m)=NaN;
                    PP_Theory_GWh{k}{l}(m)=NaN;
                    PPipeDia{k}{l}(m)=NaN;
                    OptInvP{k}{l}(m)=NaN;
                    PP_W{k}{l}(m)=NaN;
                    OptSpecCapP{k}{l}(m)=NaN;
                    OpthfP{k}{l}(m)=NaN;
                    OptDP{k}{l}(m)=NaN;
                end
            end
        end
        
        for k=1:nd
            for l = 1:numel(outlets)
                for m = 1:n_ielevations
                    DeselectedSites_unSortDP{k}{l}{m}=[];
                end
            end
        end
        
        DeselectedSites = 0;
        DeselectedSites_final = {};
        
        %%%%find_Pinlets.m END
        
        %% Determine lat lon coordinates of outlets and inlets
        disp('Calculate pipelenght and lat lon');
        
        %%%%Ppipelength.m START
        
        %% Calculate pipeplength in data window
        % Based on columns rows and elevation diff
        for  k = 1:numel(outlets)
            if isnan(ro(k)); continue; end
            for l = 1:n_ielevations
                for m = 1:numel(Pinlet_win{k}{l})
                    %PL{k}{l}(m) = Lenght_finder(acc_inlet_win{k},Zout,Pinlet_win{k}{l}(m),latOut(k),Zoutlets(k),ZPinlet{k}{l}(m));
                    PL{k}{l}(m)=Length_finderSD(acc_inlet_win{k},Poutlet_win,Pinlet_win{k}{l}(m),latOut(k),Zoutlets(k),ZPinlet{k}{l}(m),coordsys, cellsz_m);
                    
                    %Get coordinates for inlets
                    if strcmp(coordsys, 'geographic')   %get lat lon for the inlet
                        [latP{k}{l}(m), lonP{k}{l}(m)] = setltln(acc, Rw, rPin{k}{l}(m), cPin{k}{l}(m));
                    elseif strcmp(coordsys, 'planar')
                        [xP, yP] =  intrinsicToWorld(Rw,cPin{k}{l}(m), rPin{k}{l}(m));
                        [latP{k}{l}(m), lonP{k}{l}(m)] = projinv(proj,xP, yP);
                    end
                end
            end
        end
        %%%%Ppipelength.m END
        
        
        %% Identify cells that overlap w pipe to get resettlement and landuse costs
        
        for  k =  1:numel(outlets)
            if isnan(ro(k)); continue; end
            for l = 1:n_ielevations
                n_inlets=numel(Pinlet_win{k}{l});
                n_DPinlets(k)=n_DPinlets(k)+n_inlets;
                
                for m = 1:n_inlets
                    % Get cells overlapping the inlet-outlet displacement line
                    [~,~, DPcells_idx] = bresenham(cPin{k}{l}(m), rPin{k}{l}(m), co(k), ro(k), size(acc));
                    
                    %% Resettlement cost for DP
                    %Find GDP values % $/cap/year
                    [GDPr,~] = find(ISOGDP(:,1)==RDCountry_id(k));
                    if isempty(GDPr)==1; GDPpc = mean(ISOGDP(:,3)); end %if ISO value is not found, world average
                    GDPpc = ISOGDP(GDPr,3);
                    
                    % Population
                    PopDP{k}{l}(m) = sum(Pop(DPcells_idx))*(cellsz_m/1e3)^2 ;
                    MoneyDPlake = ResettleFactor * GDPpc * PopDP{k}{l}(m); %Resettlement cost
                    PopCostDP{k}{l}(m) = MoneyDPlake*(1+CommBenefitRate_DP); %PopCost = Resettlement + 10% for Benefit Sharing
                    
                    %% Land value for DP
                    if landval_set==1
                        %LandValue is $/km2. Area of 15s cell is 450m^2 x 9.09
                        %(infinite discounted yearly revenue with discount factor 10%, see perpetuity_check.m)
                        LandValuetmp = sum(LandValue(DPcells_idx))*(cellsz_m/1e3)^2 * 9;
                    else
                        LandValuetmp = 0;
                    end
                    LandValueDP{k}{l}(m) = LandValuetmp*(1+LandAcqRate_DP); %Landval = LandVal + additional fee;
                end
            end
        end
        %                LandValueDP{k}{l}(m) + PopCostDP{k}{l}(m)
        %% Running P cost model
        disp('Running pipe cost model');
        %%%%Pcostmodel.m START
        % Pipe cost model
        
        for  k = 1:n_outlets %outlets
            if isnan(ro(k)); continue; end
            for l = 1:n_ielevations %elevations
                for m = 1:numel(Pinlet_win{k}{l}) % inlets
                    
                    if isnan(lonP{k}{l}(m))==1
                        COEP{k}{l}(m)= NaN;
                        PP_energy_GWh{k}{l}(m) = NaN;
                        CostElementsP{k}{l}{m} = NaN;
                        nPipeP{k}{l}(m)=NaN;
                        continue; end
                    
                    %% Check min head and Q requirement
                    %Skip if head less than Plarge_minH or if it is in tributary and head is less than Psmall_minH
                    tmphead=ZPinlet{k}{l}(m)- Zoutlets(k);
                    if (tmphead < Plarge_minH) || (enablesmallPdesignQH && out_type(k)==1 && tmphead < Psmall_minH)
                        COEP{k}{l}(m)= NaN;
                        PPnet{k}{l}(m) = NaN;
                        CostElementsP{k}{l}{m} = NaN;
                        nPipeP{k}{l}(m)=NaN;
                        continue;
                    end
                    
                    % Skip if Qdesign <Qmin
                    if out_type(k)==1 && QDesignPinlet{k}{l}(m) < minQ_small || out_type(k)>1 &&QDesignPinlet{k}{l}(m) < minQ_large
                        COEP{k}{l}(m)= NaN;
                        PPnet{k}{l}(m) = NaN;
                        CostElementsP{k}{l}{m} = NaN;
                        nPipeP{k}{l}(m)=NaN;
                        continue;
                    end
                    
                    %% Eval P project cost
                    [COEP{k}{l}(m),PP_Theory_GWh{k}{l}(m),  PP_energy_kWh, CostElementsP{k}{l}{m}, nPipeP{k}{l}(m), OptInvP{k}{l}(m), PP_W{k}{l}(m), OptSpecCapP{k}{l}(m), OpthfP{k}{l}(m), OptDP{k}{l}(m)] = costmodel_pipesys(Zoutlets(k),ZPinlet{k}{l}(m),PL{k}{l}(m),QDesignPinlet{k}{l}(m),QDesignMeanPinlet{k}{l}(m),QDesignLFPinlet{k}{l}(m),DisTransOutlet(k),cost_constr,cost_lim,hazard_rate_RD(k),DisRoadOutlet(k), enablesmallPcost, DisSettlementOutlet(k),LandValueDP{k}{l}(m), PopCostDP{k}{l}(m), makechanges2cost);
                    
                    PP_energy_GWh{k}{l}(m) = PP_energy_kWh *1e-6; % kWh to GWh
                    
                    % Remove 0 projects
                    if COEP{k}{l}(m)==0; COEP{k}{l}(m)=NaN;end
                    if PP_energy_GWh{k}{l}(m)==0; PP_energy_GWh{k}{l}(m)=NaN;end
                    
                end
            end
        end
        
        %% Select Cheapest
        % First find lowest COE per elevation level
        for k=1:n_outlets
            
            if isnan(ro(k)); continue; end
            
            for l=1:n_ielevations
                if isempty(COEP{k}{l})==1; continue; end
                [COEPminElev{k}(l),bMinP{k}(l)] = min(COEP{k}{l}); %Lowest index per elevation level
            end
        end
        
        % Set to nan if no valid index for outlet
        for k=1:n_outlets
            
            if isnan(ro(k)); continue; end
            
            bMinP{k}(bMinP{k}==0)=NaN;
            COEPminElev{k}(COEPminElev{k}==0)=NaN;
        end
        
        % Secondly find lowest COE among elevation levels
        for k=1:n_outlets
            
            if isnan(ro(k)); continue; end
            if isempty(COEPminElev{k})==1; continue; end
            
            [~,bMinElevP(k)] = min(COEPminElev{k}); %Lowest index of all elevation levels
        end
        
        %% Transfer min info into actual data
        clear a1 a2
        for k= 1:n_outlets
            
            if isnan(ro(k)); continue; end
            if bMinElevP(k)==0; continue; end
            if isempty(bMinP{k})==1; continue; end
            
            a1(k)=bMinElevP(k); % index for min cost within each Elevation level ~ similar to l index
            a2(k)=bMinP{k}(bMinElevP(k)); % index for min cost across ni Elevation levels ~ similar to m index
            
            if a2(k)==0; continue; end
            if isnan(a2(k)); continue; end
            
            aCOEPmin(k) = COEP{k}{a1(k)}(a2(k));
            aPPnetmin(k) = PP_energy_GWh{k}{a1(k)}(a2(k));
            aPPtheorymin(k) = PP_Theory_GWh{k}{a1(k)}(a2(k));
            
            nPipeminP(k) = nPipeP{k}{a1(k)}(a2(k));
            
            
            if isempty(CostElementsP{k}{a1(k)})==1
                CostElementsPMin{k} = 0;
                dfQPmin(k)= 0;
                dfPLmin(k)= 0;
                HeadPmin(k) = 0;
                aPInletMin(k) = 0;
                aPinlet_windowMin(k) = 0;
                lat_Pin_min(k) = 0;
                lon_Pin_min(k) = 0;
                rPinMin(k) = 0;
                cPinMin(k) = 0;
                continue; end
            
            CostElementsPMin{k} = CostElementsP{k}{a1(k)}{a2(k)};
            dfQPmin(k)= QPinlet{k}{a1(k)}(a2(k));
            dfPLmin(k)= PL{k}{a1(k)}(a2(k));
            HeadPmin(k) = ZPinlet{k}{a1(k)}(a2(k)) - Zoutlets(k);
            HeadraceLPmin(k)= PL{k}{a1(k)}(a2(k)) - (ZPinlet{k}{a1(k)}(a2(k)) - Zoutlets(k));
            aPInletMin(k) = Pinlet{k}{a1(k)}(a2(k));
            aPinlet_windowMin(k) = Pinlet_win{k}{a1(k)}(a2(k));
            lat_Pin_min(k) = latP{k}{a1(k)}(a2(k));
            lon_Pin_min(k) = lonP{k}{a1(k)}(a2(k));
            rPinMin(k) = rem(aPInletMin(k)-1, nrw)+1;  %best inlet row
            cPinMin(k) = fix((aPInletMin(k)-1)/nrw)+1;  %best inlet col
            QDesignPinletMin(k) = QDesignPinlet{k}{a1(k)}(a2(k));
            QDesignLFPinletMin(k) = QDesignLFPinlet{k}{a1(k)}(a2(k));
            QDesignMeanPinletMin(k) = QDesignMeanPinlet{k}{a1(k)}(a2(k));
            ZPinletMin(k) = ZPinlet{k}{a1(k)}(a2(k));
            nPipePMin(k) = nPipeP{k}{a1(k)}(a2(k));
            OptInvPMin(k) = OptInvP{k}{a1(k)}(a2(k));
            PPMin(k) = PP_W{k}{a1(k)}(a2(k));
            OptSpecCapPMin(k) = OptSpecCapP{k}{a1(k)}(a2(k));
            accPmin(k) = accPinlet{k}{a1(k)}(a2(k));
            OpthfPMin(k) = OpthfP{k}{a1(k)}(a2(k));
            OptDPMin(k) = OptDP{k}{a1(k)}(a2(k));
            
        end
        
        %%%%Pcostmodel.m END
        
        %% Q-decreaser
        disp('Running Q decreaser for pipe systems');
        
        %%%%QPdecreaser.m START
        for k = 1:n_outlets
            %fprintf('%d DP QP decreaser #%d of %d\n',nbasin,k,n_outlets)
            
            if isnan(ro(k)); continue; end
            if isempty(aPinlet_windowMin(k)); continue; end
            if aPinlet_windowMin(k)==0; continue; end
            if isnan(aCOEPmin(k))==1; continue; end
            
            clear r c
            [r,c] = ind2sub(size(Q_inlet_win{k}),aPinlet_windowMin(k));
            
            %QQ = Qdecimater(r,c,inlet_win,Q_inlet_win{k},fdir_inlet_win{k},acc_inlet_win{k},aPinlet_windowMin(k));
            QQ = QdecimaterSD(r,c,inlet_win,Q_inlet_win{k},fdir_inlet_win{k},adir_inlet_win{k},flowdist_inlet_win{k},acc_inlet_win{k},aPinlet_windowMin(k));
            
            Q_inlet_win{k} = QQ;
        end
        
        %% Check
        % k=44;
        %
        % [rin,cin] = ind2sub(size(Q_inlet_win{k}),aPinlet_windowMin(k));
        %
        % figure(1);clf;imagesc(log(Q_inlet_win{k}));axis image; colormap(flipud(gray));
        % hold on
        % plot(inlet_win+1,inlet_win+1,'.r','markersize',20)
        % plot(cin,rin,'.b','markersize',20)
        % hold off
        
        %%%%QPdecreaser.m END
        %% Store info per Q decreaser loop
        disp('Store info per Q decreaser loop')
        
        PCOEend{ndl}=aCOEPmin; % $/kWh Diversion-Pipe COE
        PPnetend{ndl}=aPPnetmin; % GWh Diversion-Pipe energy actual potential
        PPtheoryend{ndl}=aPPtheorymin; % GWh Diversion-Pipe energy theoretical potential
        Platminend{ndl} = lat_Pin_min;
        Plonminend{ndl}= lon_Pin_min;
        rPinMinend{ndl} = rPinMin;
        cPinMinend{ndl}= cPinMin;
        aPInletminEnd{ndl}=aPInletMin;
        aPinlet_windowMinEnd{ndl}=aPinlet_windowMin;
        dfQPminEnd{ndl}=dfQPmin;
        dfPLminEnd{ndl}=dfPLmin;
        CostElementsPMinEnd{ndl}=CostElementsPMin;
        Pinlet_winend{ndl}=Pinlet_win;
        HeadPminend{ndl}=HeadPmin;
        HeadraceLPminend{ndl}=HeadraceLPmin;
        QDesignPinletMinend{ndl}=QDesignPinletMin;
        QDesignLFPinletMinend{ndl}=QDesignLFPinletMin;
        QDesignMeanPinletMinend{ndl}=QDesignMeanPinletMin;
        dfPLminend{ndl}=dfPLmin;
        ZPinletMinend{ndl} = ZPinletMin;
        nPipePMinend{ndl} = nPipePMin;
        OptInvPMinend{ndl} = OptInvPMin;
        PPMinend{ndl} = PPMin;
        OptSpecCapPMinend{ndl} = OptSpecCapPMin;
        PSurfaceLake15sMinend{ndl} = PSurfaceLake15s;
        accPminend{ndl} = accPmin;
        OpthfPminend{ndl} = OpthfPMin;
        OptDminend{ndl} = OptDPMin; %Optimal pipe diameter
        
        %%
        aDPPnetEnd4  = horzcat(PPnetend{:})'; % GWh  Pipe systems Pnet
        fprintf('>>>>>>Total P potential: %0.0f GWh\n',sum(aDPPnetEnd4(~isnan(aDPPnetEnd4))));
        %%%% PipeSystems.m END
    end
end

%% RD: CLASSIC DAM SYSTEMS with RUN-OF-RIVER CHARACTERISTICS, LOADFACTOR < 100%
if runRD==1
    fprintf('\n\n\n\nRiver dam scanner\n');
    
    %%%% riverdamscannerHR.m START
    %% Riverdam scanner to determine width and height in search radius
    %counter=0;
    %numWorkers=0;
    fprintf('Generating RD dam height-width\n')
    for  k = 1:numel(outlets)
        %parfor (k = 1:numel(outlets),numWorkers)
        % k=32;
        %fprintf('%d RD dam height-width #%d of %d\n',nbasin,k,numel(outlets))
        
        % Skip outlet if outlet is in nan
        if isnan(ro(k)); continue; end
        
        % Skip if not the index for RD
        if varyRiverLevel_HPtype && out_type(k)~=out_type_RD; continue; end
        
        % Skip if minQ not met
        if Q_RD_design(k)<minQ_large; continue; end
        
        % First load high res dem
        if strcmp(coordsys, 'geographic') % use Davids 3s tiles
            RDdemHR = hs3s_window('con',latOut(k), lonOut(k), RDwinsize_hires_dem, rootpath,continent_in);
        elseif strcmp(coordsys, 'planar')
            %% Get high-resolution window data
            [nr100,nc100]=size(Z_100m);
            %get RC from XY coordinate for 100m DEM
            [r_dis,c_dis] = worldToDiscrete(Rw_100m,xOut(k),yOut(k));
            
            %create window in high-res dem
            firstrow = max(1,r_dis-RDwinsize_hires_dem);
            lastrow  = min(nr100,r_dis+RDwinsize_hires_dem);
            firstcol = max(1,c_dis-RDwinsize_hires_dem);
            lastcol  = min(nc100,c_dis+RDwinsize_hires_dem);
            RDdemHR = Z_100m(firstrow:lastrow,firstcol:lastcol);
            
        end
        RDdemHR = single(RDdemHR);
        
        % To control for return in demHR function
        isfile_outlets(k)=0;
        if RDdemHR==0;
            isfile_outlets(k)=1;
            ro(k) = NaN;
            co(k) = NaN;
            lonOut(k)=NaN;
            latOut(k)=NaN;
            dfhdamRD{k} = 0;
            dfldamRD{k} = 0;
            ZiCenterRD(k) = 0;
            continue;
        end
        
        [nrRD,ncRD]=size(RDdemHR);
        
        % Second find lowest point in 15s highres window
        RDwin3s=2;  % Create a two cell window around center of DEM
        
        % RDwinsize_dem represents the center of the square RDdemHR, so create 2
        % cell window around it
        firstrowRD_3s = max(1,RDwinsize_hires_dem+1-RDwin3s);
        lastrowRD_3s  = min(nrRD,RDwinsize_hires_dem+1+RDwin3s);
        firstcolRD_3s = max(1,RDwinsize_hires_dem+1-RDwin3s);
        lastcolRD_3s  = min(ncRD,RDwinsize_hires_dem+1+RDwin3s);
        
        RDdemHR_3s = RDdemHR(firstrowRD_3s:lastrowRD_3s,firstcolRD_3s:lastcolRD_3s); %Subset 3s DEM tile even more
        [RDr_mindemHR,RDc_mindemHR]=find(RDdemHR_3s==min(RDdemHR_3s(:))); %Find lowest point in 2 cells window around point
        %Transformed focal point based on demHR_15s
        RDrh_width2 = RDwinsize_hires_dem+1+(RDr_mindemHR(1)-RDwin3s-1); % New coordinates in the 3s DEM window
        RDch_width2 = RDwinsize_hires_dem+1+(RDc_mindemHR(1)-RDwin3s-1); % In case multiple lowest point, select first option
        
        ZiCenterRD(k) = RDdemHR(RDrh_width2,RDch_width2); %Corrected Zfocus
        
        %Check high res DEM
        if showcheckplots
            figure(1);clf;imagesc(RDdemHR);axis image;colormap(jet);
            figure(1);hold on;
            plot(RDch_width2,RDrh_width2,'r.','markersize',20);
            hold off
            title('High res DEM window around point')
            figure(2);clf;imagesc(RDdemHR_3s);axis image;colormap(jet);
            figure(2);hold on; plot(RDc_mindemHR,RDr_mindemHR,'r.','markersize',20)
            figure(2);hold on; plot(RDc_mindemHR(1),RDr_mindemHR(1),'xk','markersize',20); hold off
            title('Zooming into immediate 2 cells window around point in high res DEM ')
        end
        
        %%Find highest point in highres window around point
        RDdemMax(k) = max(RDdemHR(:));
        
        if MiniHydro_special==1 && Q_RD_design(k) < bighydro_cutoff
            MiniFlag(k)=1;
            if (RDdemMax(k)-ZiCenterRD(k)-1) > MiniDamMax
                RDZdiff(k) = MiniDamMax;
                
            else
                RDZdiff(k) = RDdemMax(k)-ZiCenterRD(k)-1;
            end
            
        else
            MiniFlag(k)=0;
            RDZdiff(k) = RDdemMax(k)-ZiCenterRD(k)-1;
        end
        
        % Generate all possible dam heights greater than minimum dam height
        dfhdamRD{k} =RD_minH:RDZdiff(k);
        
        %% Based on dam height we determine dam width
        if ASYMRD==0
            ldamRD=0;
            
            for wRD=1:length(dfhdamRD{k})
                % Determine Ztop based on hdam and Zfocus
                clear RDZtop RDdx_w RDdy_w RDds_w RDhiground_w RDdhiground_w
                RDZtop = ZiCenterRD(k) + wRD; % get elevation of dam top
                
                % Build higround distance map
                RDdx_w = repmat(abs(-RDwinsize_hires_dem:RDwinsize_hires_dem), [2*RDwinsize_hires_dem+1, 1]);
                RDdy_w = RDdx_w';
                RDds_w = sqrt(RDdx_w.^2+RDdy_w.^2); % Distance matrix from center of window to surrounding areas
                
                % Determine minimum distance from Zfocus cell to higround
                % in the hi-res dem window
                RDhiground_w = RDdemHR>RDZtop;
                RDdhiground_w = RDds_w(RDhiground_w);
                
                ldamRD(wRD) = max(0.5,min(RDdhiground_w))*90*2; % twice for symmetric dam in m
                
                if isempty(ldamRD(wRD));ldamRD(wRD)=0;end
                
                %checkdfld
                %fprintf('lenght = %f meters.\n',ldamRD(wRD))
                %figure(1);clf;imagesc(RDhiground_w);axis image;colormap(jet);hold on
                %plot(RDch_width2,RDrh_width2,'r.','markersize',20); hold off
                
            end
            
            dfldamRD{k} = ldamRD; % dam length for given height in m
            
            
        elseif ASYMRD==1
            
            [lrwidth{k}, genwidth{k}]=AsymDamLenght(rootpath,continent_in,latOut(k),lonOut(k),RDwinsize_hires_dem,RDZdiff(k),0);
            
            dfldamRD{k} = lrwidth{k};
        end
        
    end
    %%%% riverdamscannerHR.m END
    
    %% River Dam population displacement and land loss
    disp('River dam population displacer');
    
    %%%% RDlake.m START
    % RiverDam lake to calculate population displacement
    RDupstream_arc=cell(numel(outlets),1);
    disp('Generating RD Lake Pop-Landvalue')
    for k=1:numel(outlets)
        %fprintf('%d RD Lake Pop-Landvalue #%d of %d\n',nbasin,k,numel(outlets))
        
        if isnan(ro(k)); continue; end
        % Skip if not the index for RD
        if varyRiverLevel_HPtype && out_type(k)~=out_type_RD; continue; end
        
        % Skip if minQ not met
        if Q_RD_design(k)<minQ_large; continue; end
        
        %Calculate upstream cells for each outlet
        %            RDupstream = fastfindupstream_lim(acc,fdir,drow,dcol,outIdx(k));
        RDupstream = fastfindupstream_DisSD4(acc,fdir,adir,flowdist,drow,dcol,outlets(k),0); %no plot % no dowin
        RDupstream_arc{k}=sparse(RDupstream);
        
        % define box around upstream cells to speed up
        ZupstreamRD = RDupstream==1;
        %ZupstreamRD = RDupstream==1 & Z > Z(ro(k),co(k)) & Z < (Z(ro(k),co(k))+dfhdamRD{k}(end));
        cs = find(sum(ZupstreamRD)>0); % columns
        rs = find(sum(ZupstreamRD')>0); % rows
        colwin=0;
        
        fcs = max(1,cs(1) - colwin);
        lcs = min(ncw,cs(end) + colwin);
        frs = max(1,rs(1) - colwin);
        lrs = min(nrw,rs(end) + colwin);
        Zupstream_win = Z(frs:lrs,fcs:lcs);
        RDupstream_win = RDupstream(frs:lrs,fcs:lcs);
        Popupstream_win = Pop(frs:lrs,fcs:lcs);
        LandValueupstream_win = LandValue(frs:lrs,fcs:lcs);
        
        %Initialize
        PopRDlake=0;
        MoneyRDlake=0;
        LandValueRDlakee=0;
        clear RDlake2 hd
        
        %Determine reservoir limit
        RDreslim=storage_lim* Qyravg(ro(k),co(k)) * 365*24*60*60;  % percentage of annual avg inflow volume m3/s *s
        % Loop through dam heights to evaluate res vol, area and associated
        % costs - stop loop when vol exceeds 5% of annual inflow
        for j=1:numel(dfhdamRD{k})
            hd=single(dfhdamRD{k}(j));
            
            RDlake2 = RDupstream_win & Zupstream_win < (Zoutlets(k)-RDDepth(k)+OptDH(k)); % Based on 15s DEM map, minus river depth + dam height based on 3s DEM
            %Volume
            Z_upstream = Zupstream_win(RDupstream_win);
            dz_h = Zoutlets(k)-RDDepth(k)+hd-Z_upstream;
            dz_h = max(0, dz_h);
            
            if strcmp(coordsys, 'geographic')
                %Surface
                RDlakeSurface{k}(j) = sum(RDlake2(:)==1)* (cellsz_m*(cellsz_m*cosd(latOut(k)))); % m2
                %Volume
                RDVolumeLake{k}(j)  = sum(dz_h * cellsz_m * (cellsz_m*cosd(latOut(k)))); % m3
            elseif strcmp(coordsys, 'planar')
                %Surface
                RDlakeSurface{k}(j) = sum(RDlake2(:)==1)* cellsz_m^2; % m2
                %Volume
                RDVolumeLake{k}(j)  = sum(dz_h * cellsz_m^2); % m3
            end
            % Find indices of lake cellsl
            RDlakeidx = find(RDlake2);
            
            %% Break loop and shorten dam height options if reservoir exceeds storage limit
            if enablereslimit && RDVolumeLake{k}(j)>RDreslim
                dfhdamRD{k}=dfhdamRD{k}(1:j-1);
                dfldamRD{k}=dfldamRD{k}(1:j-1);
                RDlakeSurface{k} =RDlakeSurface{k}(1:j-1);
                RDVolumeLake{k}=RDVolumeLake{k}(1:j-1);
                %disp({'Break at ' j})
                break
            end
            
            %% Find GDP and resettlement cost
            %GDPSwiss = 55000; % $/cap/year
            %GDPUSA = 54000; % $/cap/year
            [GDPr,~] = find(ISOGDP(:,1)==RDCountry_id(k));
            if isempty(GDPr)==1; GDPpc = mean(ISOGDP(:,3)); end %if ISO value is not found, world average
            GDPpc = ISOGDP(GDPr,3);
            
            PopRDlake(j) = sum(Popupstream_win(RDlakeidx))*(cellsz_m/1e3)^2 ; %inhabitants /km2 * cell area in km
            MoneyRDlake(j) = ResettleFactor * GDPpc * PopRDlake(j); %Resettlement cost
            
            %% Find land value costs
            if landval_set==1
                %LandValue is $/km2. Area of 15s cell is 450m^2 x 9.09
                %(infinite discounted yearly revenue with discount factor 10%, see perpetuity_check.m)
                LandValueRDlakee(j) = sum(LandValueupstream_win(RDlakeidx))*(cellsz_m/1e3)^2 * 9.09;
            else
                LandValueRDlakee(j) = 0;
            end
            
        end
        
        % Save resettlement and land value costs
        PopDisplaced{k}=PopRDlake;
        %BenefitCost=10/100*MoneyRDlake;
        PopCost{k}=MoneyRDlake*(1+CommBenefitRate_RP); %PopCost = Resettlement + 10% for Benefit Sharing
        %RDLake_map{k}=RDlake2;
        LandValueRDlake{k}=LandValueRDlakee*(1+LandAcqRate_RP); %Landval = LandVal + additional fee;
        
        %% Check dam heights vs lengths and reservoir volumes
        if showcheckplots
            figure;
            subplot(1,2,1)
            plot(dfhdamRD{k},dfldamRD{k})
            xlabel("Dam heights (m)");   ylabel("Dam widths (m)")
            grid on
            subplot(1,2,2)
            plot(dfhdamRD{k},RDVolumeLake{k})
            hold on
            yline(RDreslim,'r:')
            xlabel("Dam heights (m)");   ylabel("Reservoir volume (m^3)")
            set(gca, 'YScale', 'log')
            grid on
        end
    end
    
    %% Check lake
    % figure(1);clf;
    % img1= truecolorsc(Zupstream_win,flipud(gray));
    % % img2 = burnmask(img1, ~RDLake_map{14}{40});
    % img2 = burnmask(img1, ~RDlake2);
    % %             img2 = burnmask(img1, ~RDLake_map{k}{200});
    % image(img2); axis image;
    % title('Discharge','FontSize',14);
    % %     axis([75 125 75 125]);
    % set(gca,'xtick',[],'ytick',[])
    % hold on;
    % % plot(co(23),ro(23),'r.','markersize',15);
    % hold off
    %
    %
    %[nr nc] = size(Zupstream_win)
    % for r=1:nr
    %     for c=1:nc
    %         if RDupstream_win(r,c)==1
    %             Zup(r,c) = Zupstream_win(r,c);
    %         else
    %             Zup(r,c) = NaN;
    %         end
    %     end
    % end
    % figure(3);clf;imagesc(Zup); axis image
    %
    % %%
    % for r=1:nr
    %     for c=1:nc
    %         Zupm(r,c) = 710 - Zup(r,c);
    %     end
    % end
    % figure(4);clf;imagesc(Zupm); axis image
    %
    % %%
    % Zupmm = max(0, Zupm);
    % figure(5);clf;imagesc(Zupmm); axis image
    
    %%%% RDlake.m END
    
    %% River Dam scanner
    disp('River dam Cost model');
    
    %%%%CostModel_RiverDam.m START
    %% Costmodel_RiverDam
    % Calculate cost and potential of river dam
    %     COETotRD = 0;
    %     RDPnet=0;
    clear DH
    for  k = 1:n_outlets
        % Skip outlet if outlet is in nan
        if isnan(ro(k))
            COETotRD(k)=0;
            RDPnet(k)=0;
            RDtheory_GWh(k)=0;
            OptDH(k)=0;
            OptDL(k)=0;
            RDCostElements{k}=0;
            continue;
        end % Skip empty arrays
        
        % Skip if not the index for RD
        if varyRiverLevel_HPtype && out_type(k)~=out_type_RD
            COETotRD(k)=0;
            RDPnet(k)=0;
            RDtheory_GWh(k)=0;
            OptDH(k)=0;
            OptDL(k)=0;
            RDCostElements{k}=0;
            continue; end
        
        % Skip if minQ not met
        if Q_RD_design(k)<minQ_large
            COETotRD(k)=0;
            RDPnet(k)=0;
            RDtheory_GWh(k)=0;
            OptDH(k)=0;
            OptDL(k)=0;
            RDCostElements{k}=0;
            continue;
        end % Skip empty arrays
        
        if numel(dfldamRD{k})==0
            COETotRD(k)=0;
            RDPnet(k)=0;
            RDtheory_GWh(k)=0;
            OptDH(k)=0;
            OptDL(k)=0;
            RDCostElements{k}=0;
            continue;
        end % Skip empty arrays
        
        if numel(dfhdamRD{k})==0
            COETotRD(k)=0;
            RDPnet(k)=0;
            RDtheory_GWh(k)=0;
            OptDH(k)=0;
            OptDL(k)=0;
            RDCostElements{k}=0;
            continue;
        end % Skip empty arrays
        
        DH= single(dfhdamRD{k});
        
        if MiniHydro_special==1 && Q_RD_design(k) < bighydro_cutoff
            [COETotRD(k), RDPnet(k), OptP(k), OptDH(k), OptDL(k), RDCostElements{k}]=MiniCostmodel(dfldamRD{k},DH,PopCost{k},Q_RD_design(k),Q_RD_design_LF(k),DisTransOutlet(k),LandValueRDlake{k},RDDepth(k),nbasin,k,cost_constr,cost_lim);
        else
            [COETotRD(k), RDtheory_GWh(k), RDPnet(k), OptP(k), OptDH(k), OptDL(k), RDCostElements{k}, OptInv(k), RDP(k), OptPop(k), OptLV(k), OptSpecCap(k)]=RDcostmodel(dfldamRD{k},DH,PopCost{k},Q_RD_design(k),Q_RD_design_LF(k),DisTransOutlet(k),LandValueRDlake{k},RDDepth(k),nbasin,k,cost_constr,cost_lim,hazard_rate_RD(k),DisRoadOutlet(k), makechanges2cost);
        end
        idx_OptDH=find( DH==OptDH(k));
        PopDisplacedOpt(k) = PopDisplaced{k}(idx_OptDH);
    end
    
    
    % Remove zeros
    for k =1:n_outlets
        if COETotRD(k)==0; COETotRD(k)=NaN; end
        if RDPnet(k)==0; RDPnet(k)=NaN;end
    end
    
    COETotRD = COETotRD';
    RDPnet = RDPnet'; % kWh
    RDPnet = RDPnet*1e-6; % convert to GWh
    %%%%CostModel_RiverDam.m END
    
    %% River dam selector
    fprintf('\n\n\n Dam selector\n');
    if runMini==0
        %%%% Dam_selector.m START
        %% Dam selector
        % Routine to de-selects dams according to a certain priority
        
        if sum(isnan(COETotRD))==numel(ro)
            RDlake_Opt=0;
            CheapestDam=0;
            DeselectedSites=0;
            DeselectedSites_unSort=0;
            COETotRDs=0;
            RDPnets=0;
        end
        
        %% Collecting all dam locations DP, P and RD
        if runDP==1
            %Zeros are NaN
            aInletminEnd{1}(aInletminEnd{1}(:)==0)=NaN;
            aInletminEnd{2}(aInletminEnd{2}(:)==0)=NaN;
            outlets(outlets(:)==0)=NaN;
            
            DPlocs = horzcat(aInletminEnd{:})';
            Plocs  = horzcat(aPInletminEnd{:})';
            RDlocs = outlets';
            
            Damlocs = [DPlocs;Plocs;RDlocs]; %In terms of basin indices
            
            % Vector to keep track of type of systems
            SysID  = zeros((numel(RDlocs)+numel(DPlocs)+numel(Plocs)),1);
            SysID(1:numel(DPlocs))=1;                                      %DP systems
            SysID((numel(DPlocs)+1):(numel(DPlocs)+numel(Plocs)))=2;       %P systems
            SysID((numel(DPlocs)+numel(Plocs)+1):end)=3;                   %RD systems
            
            % Rows and columns
            ros = [horzcat(rinMinend{:})'; horzcat(rPinMinend{:})'; ro];
            cos = [horzcat(cinMinend{:})'; horzcat(cPinMinend{:})'; co];
            
            lats = [horzcat(latminend{:})'; horzcat(Platminend{:})'; latOut];
            lons = [horzcat(latminend{:})'; horzcat(Plonminend{:})'; lonOut];
            
            % %%
            % figure(1);clf;imagesc(log(Q));colormap(flipud(gray));axis image;
            % hold on
            % plot(co(find(~isnan(COETotRD))),ro(find(~isnan(COETotRD))),'.r','markersize',20)
            % plot(cinMinend{1}(find(~isnan(COEend{1}))),rinMinend{1}(find(~isnan(COEend{1}))),'.b','markersize',20)
            % plot(cPinMinend{1}(find(~isnan(PCOEend{1}))),rPinMinend{1}(find(~isnan(PCOEend{1}))),'.g','markersize',20)
            % hold off
            
        else
            %Zeros are NaN
            outlets(outlets(:)==0)=NaN;
            
            Plocs  = horzcat(aPInletminEnd{:})';
            RDlocs = outlets;
            
            Damlocs = [Plocs;RDlocs]; %In terms of basin indices
            Grandlocs = GrandIdx; %Existing locations
            
            % Vector to keep track of type of systems
            SysID  = zeros((numel(RDlocs)+numel(Plocs)),1);
            SysID((1:numel(Plocs)))=1;       %P systems
            SysID((numel(Plocs)+1):end)=2;   %RD systems
            
            % Rows and columns
            ros = [horzcat(rPinMinend{:})'; ro];
            cos = [horzcat(cPinMinend{:})'; co];
            
            lats = [horzcat(Platminend{:})'; latOut'];
            lons = [horzcat(Plonminend{:})'; lonOut'];
            
        end
        
        %% Show
        % [rRD,cRD]=ind2sub(size(acc),RDlocs);
        % % [rDP,cDP]=ind2sub(size(acc),DPlocs);
        % [rP,cP]=ind2sub(size(acc),Plocs);
        %
        % figure(1);clf;imagesc(log(Q));axis image;colormap(flipud(gray));
        % hold on
        % % plot(cRD,rRD,'r.','markersize',10)
        % % plot(cDP,rDP,'b.','markersize',20)
        % % plot(cP,rP,'g.','markersize',10)
        % plot(coss,ross,'g.','markersize',20)
        % plot(c_damst,r_damst,'b.','markersize',20)
        % hold off
        %% Collecting COEs and Pnets
        if runDP==1
            % COEs
            %             COEendDP = horzcat(COEend{:});
            %             COEendDP(COEendDP(:)==0)=NaN;
            %
            %             COEendP = horzcat(PCOEend{:});
            %             COEendP(COEendP(:)==0)=NaN;
            %
            %             COEAll = [COEendDP';COEendP';COETotRD];
            %
            %             % Pnets
            %             PnetendDP = horzcat(DPPnetend{:});
            %             PnetendDP(PnetendDP(:)==0)=NaN;
            %
            %             PnetendP = horzcat(PPnetend{:});
            %             PnetendP(PnetendP(:)==0)=NaN;
            %
            %             PnetAll = [PnetendDP';PnetendP';RDPnet];
            
        else
            % COEs
            COEendP = horzcat(PCOEend{:});
            COEendP(COEendP(:)==0)=NaN;
            
            COEAll = [COEendP';COETotRD];
            
            % Pnets
            PnetendP = horzcat(PPnetend{:});
            PnetendP(PnetendP(:)==0)=NaN;
            
            PnetAll = [PnetendP';RDPnet];
            
        end
        
        %% First, recalculate lake based on optimal dam height and check which lake floods which dams
        %% DP-systems --- NOT RUN RIGHT NOW as runDP=0
        if runDP==1
            if sum(isnan(COEendDP))~=numel(COEendDP)
                for ndd=1:2
                    for k=1:n_outlets
                        fprintf('Recalculating DP lakes outlet #%d of %d\n',k,n_outlets);
                        
                        clear DPupstream a OutletDeselect
                        if isnan(ro(k));
                            DeselectedSites_unSortDP{ndd}{k}=0;
                            continue; end
                        if isnan(aInletminEnd{ndd}(k))==1;
                            DeselectedSites_unSortDP{ndd}{k}=0;
                            continue; end
                        
                        %DPupstream = fastfindupstream_lim(acc,fdir,drow,dcol,aInletminEnd{ndd}(k));
                        DPupstream = fastfindupstream_DisSD4(acc,fdir,adir,flowdist,drow,dcol,aInletminEnd{ndd}(k),0); %no plot % no dowin
                        
                        DPlake_Opt = DPupstream & Z < (Z(aInletminEnd{ndd}(k))-DPDepth{ndd}(k)+ahdammin{ndd}(k));
                        
                        LakeIdx = find(DPlake_Opt);
                        a = ismember(Damlocs,LakeIdx);
                        OutletDeselect = find(a);
                        kc=k;
                        if ndd==2; kc=k+n_outlets;end % Making sure it doesn't select itself in the second round
                        DeselectedSites_unSortDP{ndd}{k} = OutletDeselect(OutletDeselect~=kc); % Which dams are flooded (unsorted)
                        
                        clear a
                        a = ismember(Grandlocs,LakeIdx);
                        GrandOutletDeselectDO = find(a);
                        
                    end
                end
            end
        end
        
        %% Show stuation before deselection
        
        % for i=1:numel(DPlake_Opt{1})
        %     a(i)=isempty(DPlake_Opt{1}{i});
        % end
        % b=find(~a);
        %
        % %
        % figure(1);clf;
        % img{1}= truecolorsc(DPlake_Opt{1}{b(1)},flipud(gray));
        % %
        % for i=1:numel(b)
        %     if i==1; continue; end
        %         img{i} = burnmask(img{i-1}, ~DPlake_Opt{1}{b(i)});
        % end
        % image(img{end}); axis image;
        % %
        % hold on
        % plot(cinMinend{1}(:),rinMinend{1}(:),'.r','markersize',15)
        % plot(cPinMinend{1}(:),rPinMinend{1}(:),'.g','markersize',15)
        % plot(cinMinend{2}(:),rinMinend{2}(:),'.r','markersize',15)
        % plot(cPinMinend{2}(:),rPinMinend{2}(:),'.g','markersize',15)
        % plot(co,ro,'.b','markersize',15)
        % hold off
        
        
        %% RD-systems
        fprintf('# of valid RD-systems: %d of %d\n',sum(isnan(COETotRD)),n_outlets);
        if sum(isnan(COETotRD))~=numel(ro)
            for k=1:n_outlets
                %fprintf('Recalculating lakes outlet #%d of %d\n',k,n_outlets);
                
                clear a OutletDeselect Z_upstream RDupstream dz_h RDlake_Opt
                % Skip if not the index for RD
                if varyRiverLevel_HPtype && out_type(k)~=out_type_RD
                    
                    DeselectedSites_unSortRD{k}=0;
                    DeselectedGrandSites_unSortRD{k}=0;
                    continue; end
                
                if isnan(ro(k))
                    DeselectedSites_unSortRD{k}=0;
                    DeselectedGrandSites_unSortRD{k}=0;
                    continue; end
                if isnan(COETotRD(k))==1
                    DeselectedSites_unSortRD{k}=0;
                    DeselectedGrandSites_unSortRD{k}=0;
                    continue; end
                
                
                %RDupstream = fastfindupstream_lim(acc,fdir,drow,dcol,outlets(k));
                RDupstream = RDupstream_arc{k};
                %RDupstream = fastfindupstream_DisSD5(acc,fdir,adir,flowdist,drow,dcol,outlets(k),0); %no plot % no dowin
                
                %RDlake2 = RDupstream_win & Zupstream_win < (Zoutlets(k)-RDDepth(k)+OptDH(k)); % Based on 15s DEM map, minus river depth + dam height based on 3s DEM
                
                % Based on 15s DEM map, minus river depth + dam height based on 3s DEM
                RDlake_Opt = RDupstream & Z < (Zoutlets(k)-RDDepth(k)+OptDH(k));
                RDlake_Opt_arc{k}= RDlake_Opt;
                
                %For potential dams
                LakeIdx = find(RDlake_Opt);
                a = ismember(Damlocs,LakeIdx);
                OutletDeselect = find(a);
                if runDP==1
                    DeselectedSites_unSortRD{k} = OutletDeselect(OutletDeselect~=(numel(DPlocs)+numel(Plocs)+k)); % Which dams are flooded (unsorted)
                else
                    DeselectedSites_unSortRD{k} = OutletDeselect(OutletDeselect~=(numel(Plocs)+k)); % Which dams are flooded (unsorted)
                end
                
                %For existing dams
                clear a
                a = ismember(Grandlocs,LakeIdx);
                GrandOutletDeselectRD = find(a); %
                DeselectedGrandSites_unSortRD{k} = GrandOutletDeselectRD; % Which Grand dams are flooded (unsorted)
                
                %% Lake volume calculator m3
                Z_upstream = Z(RDupstream);
                dz_h = Zoutlets(k)-RDDepth(k)+OptDH(k)-Z_upstream;
                dz_h = max(0, dz_h);
                
                if strcmp(coordsys, 'geographic')
                    RDVolumeLake15s(k)  = sum(dz_h * cellsz_m * (cellsz_m*cosd(latOut(k)))); %Volume lake m3
                    RDSurfaceLake15s(k) = max(1,sum(RDlake_Opt(:))) * cellsz_m * (cellsz_m*cosd(latOut(k))); %Surface Reservoir m2
                    
                elseif strcmp(coordsys, 'planar')
                    RDVolumeLake15s(k)  = sum(dz_h * cellsz_m * cellsz_m); %Volume lake m3
                    RDSurfaceLake15s(k) = max(1,sum(RDlake_Opt(:))) * cellsz_m^2; %Surface Reservoir m2
                end
            end
            
        end
        
        %% Collecting all deselected DamLocs
        %clear DeselectedSites CheapestDam DeselectedSites_final
        
        if runDP==1
            DeselectedSites_unSortDPP = horzcat(DeselectedSites_unSortDP{:});
            
            for k = 1:numel(DeselectedSites_unSortDPP); DeselectedSites_unSortP{k}=[]; end %Creating empty P array (floods nothing)
            
            DeselectedSites_unSort = [DeselectedSites_unSortDPP,DeselectedSites_unSortP,DeselectedSites_unSortRD];
            
        else
            for k = 1:(nd*numel(DeselectedSites_unSortRD)); DeselectedSites_unSortP{k}=[]; end %Creating empty P array (floods nothing)
            
            DeselectedSites_unSort = [DeselectedSites_unSortP,DeselectedSites_unSortRD];
            DeselectedGrandSites_unSort = [DeselectedSites_unSortP,DeselectedGrandSites_unSortRD];
            
            % Lake surfaces
            LakeSurfacesAll = [horzcat(PSurfaceLake15sMinend{:})'; RDSurfaceLake15s'];
            
            % Flow accumulations
            accAll = [horzcat(accPminend{:})'; Accoutlets'];
            
            % Basin IDs
            BIDall = zeros(size(accAll));
            BIDall(:) = nbasin;
            
            % Continent IDs
            CIDall = zeros(size(accAll));
            CIDall(:) = CID;
        end
        
        %% Show stuation before deselection
        % figure(1);clf;
        %
        % img{1}= truecolorsc(RDlake_Opt{1},flipud(gray));
        %
        % for i=1:numel(RDlake_Opt)
        %     if i==1; continue; end
        %
        %     img{i} = burnmask(img{i-1}, ~RDlake_Opt{i});
        % end
        % image(img{end}); axis image;
        %
        % hold on
        % plot(co(1:i),ro(1:i),'.r','markersize',15)
        % hold off
        
        %% Second, sorting through indexing
        
        if sum(isnan(COETotRD))~=numel(ro)
            if exist('DeselectedSites_unSort')
                
                outletsRel=find(~isnan(COEAll)); %Find relevant outlet idx
                
                %Cost Priority
                minCOE = min(COEAll);
                COEIndexing = 1./(COEAll/minCOE); %Lowest COE gets 1
                if sum(isnan(COEIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in COE deselector');return;end
                if sum(isinf(COEIndexing(outletsRel))) > 0; returnID=1;disp('Infs in COE deselector');return;end
                
                %Power priority
                maxPnet = max(PnetAll);
                PnetIndexing = PnetAll/maxPnet; %Highest Pnet gets 1
                if sum(isnan(PnetIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in Pnet deselector');return;end
                if sum(isinf(PnetIndexing(outletsRel))) > 0; returnID=1;disp('Infs in Pnet deselector');return;end
                
                %Cost per power priority
                COEperPnet = COEAll./PnetAll;
                minCOEperPnet = min(COEperPnet);
                COEperPnetIndexing = 1./(COEperPnet/minCOEperPnet); %Lowest COEperPnet gets 1
                if sum(isnan(COEperPnetIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in COEPnet deselector');return;end
                if sum(isinf(COEperPnetIndexing(outletsRel))) > 0; returnID=1;disp('Infs in COEPnet deselector');return;end
                
                %Lake surface priority
                MaxSurfaceLake = max(LakeSurfacesAll);
                LakeIndexing = 1./(LakeSurfacesAll/MaxSurfaceLake); %Lowest surface gets 1
                if sum(isnan(LakeIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in Lake deselector');return;end
                if sum(isinf(LakeIndexing(outletsRel))) > 0; returnID=1;disp('Infs in Lake deselector');return;end
                
                %Flow accumulation priority
                MaxaccAll = max(accAll);
                AccIndexing = 1./(accAll/MaxaccAll); %Lowest flow acc gets 1
                if sum(isnan(AccIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in acc deselector');return;end
                if sum(isinf(AccIndexing(outletsRel))) > 0; returnID=1;disp('Infs in acc deselector');return;end
                
                WeightedCOEPnet = (COEIndexing*betaC)+(PnetIndexing*betaP)+(COEperPnetIndexing*betaCP)+(LakeIndexing*betaL)+(AccIndexing*betaA);
                
                %Sorting
                [~, sortIdx] = sort(WeightedCOEPnet(outletsRel),'descend');
                
                CheapestDam = outletsRel(sortIdx);
                DeselectedSites = DeselectedSites_unSort(outletsRel(sortIdx)); %Sorted flooded dams
                
            else
                RDlake_Opt=0;
                CheapestDam=0;
                DeselectedSites=0;
                COEAll=0;
                PnetAll=0;
            end
        end
        
        %% Third, deselection starting with cheapest
        if iscell(DeselectedSites)
            CheapestDam2=CheapestDam;
            for i=1:numel(CheapestDam2)
                if sum(ismember(CheapestDam2(i+1:end), DeselectedSites{i}))>=1
                    idxdeselect{i} = find(ismember(CheapestDam2(i+1:end), DeselectedSites{i})) + i;
                    CheapestDam2(idxdeselect{i})=1;
                    continue; end
            end
        end
        
        %% Fourth, deselection starting with most expensive
        if iscell(DeselectedSites)
            CheapestDam3=CheapestDam2;
            for i=0:numel(CheapestDam3)-1
                if i==numel(CheapestDam3)-1; continue; end
                if sum(ismember(CheapestDam3(:), DeselectedSites{numel(CheapestDam3)-i})) >=1
                    CheapestDam3(numel(CheapestDam3)-i)=1;
                end
            end
        end
        
        %% Final check, none of the dams should be flooded
        if iscell(DeselectedSites)
            for i=1:numel(CheapestDam)
                if CheapestDam3(i)==1; continue; end
                DeselectedSites_final{i} = DeselectedSites{i};
            end
            
            for i=1:numel(DeselectedSites_final)
                if isrow(DeselectedSites_final{i})==0; DeselectedSites_final{i}=DeselectedSites_final{i}'; end
            end
            
            
            if exist('DeselectedSites_final')
                DeselectedSites_Final_vector = horzcat(DeselectedSites_final{:})';
                FloodCheck = sum(ismember(CheapestDam3(:),DeselectedSites_Final_vector(:)));
                fprintf('Flood check = %0.0f\n',FloodCheck);
            end
        end
        
        %% Check which index gives problem
        % for i=1:numel(CheapestDam3)
        %     x=sum(find(DeselectedSites_final{i}==899));
        %     if x>=1
        %         i
        %     end
        % end
        
        %% New COETotRD and RDPnet
        COEAlls = zeros(size(COEAll));
        COEAlls(COEAlls==0)=NaN;
        PnetAlls = zeros(size(PnetAll));
        PnetAlls(PnetAlls==0)=NaN;
        ross = zeros(size(ros));
        ross(ross==0)=NaN;
        coss = zeros(size(cos));
        coss(coss==0)=NaN;
        SysIDs = zeros(size(SysID));
        SysIDs(SysIDs==0)=NaN;
        latss = zeros(size(lats));
        latss(latss==0)=NaN;
        lonss = zeros(size(lons));
        lonss(lonss==0)=NaN;
        LakeSurfacesAlls = zeros(size(LakeSurfacesAll));
        LakeSurfacesAlls(LakeSurfacesAlls==0)=NaN;
        accAlls = zeros(size(accAll));
        accAlls(accAlls==0)=NaN;
        BIDs = zeros(size(BIDall));
        BIDs(BIDs==0)=nbasin;
        CIDs = zeros(size(CIDall));
        CIDs(CIDs==0)=CID;
        
        if iscell(DeselectedSites)
            for i=1:numel(CheapestDam3)
                if CheapestDam3(i)>1
                    COEAlls(CheapestDam3(i))=COEAll(CheapestDam3(i));
                    PnetAlls(CheapestDam3(i))=PnetAll(CheapestDam3(i));
                    
                    SysIDs(CheapestDam3(i))=SysID(CheapestDam3(i));
                    ross(CheapestDam3(i))=ros(CheapestDam3(i));
                    coss(CheapestDam3(i))=cos(CheapestDam3(i));
                    
                    latss(CheapestDam3(i))=lats(CheapestDam3(i));
                    lonss(CheapestDam3(i))=lons(CheapestDam3(i));
                    
                    LakeSurfacesAlls(CheapestDam3(i))=LakeSurfacesAll(CheapestDam3(i));
                    
                    accAlls(CheapestDam3(i))=accAll(CheapestDam3(i));
                    BIDs(CheapestDam3(i))=BIDall(CheapestDam3(i));
                    CIDs(CheapestDam3(i))=CIDall(CheapestDam3(i));
                end
            end
            
            fprintf('Total potential before deselection: %0.0f GWh\n',sum(PnetAll(~isnan(PnetAll))));
            fprintf('Total potential after deselection: %0.0f GWh\n',sum(PnetAlls(~isnan(PnetAlls))));
            
            %Deselecting dams that flood existing Grand dams
            for i=1:numel(DeselectedGrandSites_unSort)
                if isempty(DeselectedGrandSites_unSort{i})==1; continue; end
                for j=1:numel(DeselectedGrandSites_unSort{i})
                    %fprintf('%d %d. %0.2f\n',i,DeselectedGrandSites_unSort{i}(j),COEAlls(i))
                    COEAlls(i)=NaN;
                    PnetAlls(i)=NaN;
                    SysIDs(i)=NaN;
                    ross(i)=NaN;
                    coss(i)=NaN;
                    latss(i)=NaN;
                    lonss(i)=NaN;
                    LakeSurfacesAlls(i)=NaN;
                    accAlls(i)=NaN;
                    BIDs(i)=NaN;
                    CIDs(i)=NaN;
                end
            end
            fprintf('Total potential after deselecting Grand: %0.0f GWh\n',sum(PnetAlls(~isnan(PnetAlls))));
            
            DPIDs = numel(find(SysIDs==1));
            PIDs = numel(find(SysIDs==2));
            RDIDs = numel(find(SysIDs==3));
        end
        
        %% Show dams and lake
        %         figure(2)
        %         clear img
        %         img{1}= truecolorsc(RDlake_Opt_arc{CheapestDam3(1)},flipud(gray));
        %         for i=1:numel(CheapestDam3)
        %             if i==1; continue; end
        %             if CheapestDam3(i)==1;
        %                 img{i} = img{i-1};
        %                 continue;
        %             end
        %
        %             img{i} = burnmask(img{i-1}, ~RDlake_Opt_arc{CheapestDam3(i)});
        %         end
        %         %
        %         image(img{end}); axis image;
        %         title('Second selection')
        %         hold on
        %         plot(co(CheapestDam3),ro(CheapestDam3),'.r','markersize',15)
        %         hold off
        
        %         %% Show top-tot-teen graph
        %         % Dams can still overlap when they're in different streams
        %         indices_select = CheapestDam(find(~isnan(ro(CheapestDam3))));
        %         figure(5); clf;
        %         hold on;
        %         for i=1:numel(indices_select)
        %             zbot = Z(ro(indices_select(i)),co(indices_select(i)));
        %             ztop = zbot + OptDH(indices_select(i));
        %             carea = acc(ro(indices_select(i)),co(indices_select(i)));
        %             x = [log(carea) log(carea)];
        %             y = [zbot ztop];
        %             plot(x,y, 'linewidth',2,'color',rand(1,3));
        %         end
        %         hold off;
        %         xlabel('log(area)');
        %         ylabel('Elevation');
        %
        %         [lo,la] = setltln(acc, Rw, r_damst(10), c_damst(10)); %Coordinates outlets
        %
        %         %%%% Dam_selector.m END
    end
    
    % catch error
    %if returnID; continue; end
    
end

%% Report
if runMini==0 && runP==0
    fprintf('Tot potential: %0.0f GWh\n',sum(PnetAlls(~isnan(PnetAlls))));
elseif runP==1
    fprintf('\n\n\n\nTot PnetAlls potential selected: %0.0f GWh\n',sum(PnetAlls(~isnan(PnetAlls))));
    fprintf('Total P potential: %0.0f GWh\n',sum(aDPPnetEnd4(~isnan(aDPPnetEnd4))));
    fprintf('Tot RD potential initial shortlist: %0.0f GWh\n',sum(RDPnet(~isnan(RDPnet))));
    fprintf('Tot RD + P potential initial shortlist: %0.0f GWh\n',sum(RDPnet(~isnan(RDPnet)))+sum(aDPPnetEnd4(~isnan(aDPPnetEnd4))));
else
    fprintf('Tot RD potential: %0.0f GWh\n',sum(RDPnet(~isnan(RDPnet))));
end

%     catch ME %Error handling
%
%        fprintf('\nError in basin %d:\n\n%s\n\n',nbasin, ME.message);
%
%        errormessage = sprintf('%s',ME.message);
%
%        errpath = fullfile(root, sprintf('output\\%s\\%s',version, continent_in));
%        if ~isdir(errpath)
%            mkdir(errpath);
%        end
%
%        txtfile = fullfile(root, sprintf('output\\%s\\%s\\Error_%d.txt',version, continent_in, nbasin));
%        fileID = fopen(txtfile,'w');
%        fprintf(fileID,'%s',errormessage);
%        fclose(fileID);

%   end

%% SD summary
fprintf('\n\nFor %s\n', runname)
fprintf('# of initial outlets evaluated: %0.0f out of %0.0f\n',n_outlets-n_outlets_1stdeselect,n_outlets)
fprintf('%s Technical Potential: %0.2f GWh\n',scenario,nansum(PnetAlls))
fprintf('# of selected projects: %0.0f out of %0.0f\n',sum(~isnan(PnetAlls)),length(PnetAlls))
fprintf('# of selected River type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==2),sum(PnetAlls(SysIDs==2)))
fprintf('# of selected Diversion type projects: %0.0f with %0.2f GWh\n\n',sum(SysIDs==1),sum(PnetAlls(SysIDs==1)))

fprintf('%s Financial Potential: %0.2f GWh\n',scenario, nansum(PnetAlls(COEAlls<=0.1)))
fprintf('# of selected projects: %0.0f out of %0.0f\n',sum(~isnan(PnetAlls(COEAlls<=0.1))),length(PnetAlls))
fprintf('# of selected River type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==2 & COEAlls<=0.1),sum(PnetAlls(SysIDs==2 & COEAlls<=0.1)))
fprintf('# of selected Diversion type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==1 & COEAlls<=0.1),sum(PnetAlls(SysIDs==1 & COEAlls<=0.1)))

%% Display map
[rch,cch]=ind2sub(size(channel_main_trib), find(channel_main_trib>0));

f1=figure(100*labindex);clf;
inbasin=single(~outside);
inbasin(outside)=nan;
mksz=10;
imagescnan(inbasin);axis image;%colormap(gray)
hold on
legend('-DynamicLegend','location','Best')
plot(cch,rch,'.b','markersize',2,'DisplayName',"Channel") %plot channel
if ExistingDams_constraint; plot(c_damst,r_damst,'+','markersize',6,'DisplayName',"Existing Dams");end
plot(co_arc,ro_arc,'.','markersize',mksz-2,'Color',0.7*[1 1 1],'DisplayName','Initial sites');
plot(cos,ros,'o','markersize',mksz-2,'DisplayName',"Sites where cost was evaluated") %,'Color',0.7*[1 1 1]
plot(coss(SysIDs==1),ross(SysIDs==1),'.g','markersize',mksz,'DisplayName',"Optimal Diversion Projects") %inlet
plot(coss(SysIDs==2),ross(SysIDs==2),'.r','markersize',mksz,'DisplayName',"Optimal River Type Projects") %outlet

%plot(co,ro,'o','markersize',mksz,'DisplayName',"All Outlets Evaluated") %,'Color',0.7*[1 1 1]
plot(co_arc(isnan(co)),ro_arc(isnan(co)),'o','markersize',mksz-2,'DisplayName','Deselected sites');
grid on
title(sprintf('Evaluated and optimal projects in %s\nTechnical Potential: %0.2f TWh \nFinancial Potential: %0.2f TWh',runname,nansum(PnetAlls)/1000,nansum(PnetAlls(COEAlls<=0.1)/1000)),'Interpreter','none')

% if show 10 biggest dams
PnetAllsn = PnetAlls(~isnan(PnetAlls));
rossn = ross(~isnan(PnetAlls));
cossn = coss(~isnan(PnetAlls));
[~,idx] = sort(PnetAllsn,'descend');
plot(cossn(idx(1:5)),rossn(idx(1:5)),'or','markersize',mksz+3,'LineWidth',2,'DisplayName','Largest 5 dam sites');

% if GCS then apply ylim to remove blank space
ylim([0 2000])

%% Save and sum
if savedata
    disp('Saving energy data')
    %     matpath = fullfile(root, sprintf('output\\%s', continent_in));
    %     if ~isdir(matpath)
    %         mkdir(matpath);
    %     end
    saveas(f1, strrep(matfile,'.mat','.fig'))
    %% Save summary tables in Excel
    PID=transpose(1:numel(lats)); %create ID for all feasible locations
    outdata_all=table(PID,lats,lons,SysID,COEAll,PnetAll); % selection after DP and R searches
    outdata_ss=table(PID,latss,lonss,SysIDs,COEAlls,PnetAlls); % selection after overlap check
    
    writetable(outdata_ss,xlsfile,'Sheet','outdata_ss')
    writetable(outdata_all,xlsfile,'Sheet','outdata_all')
    %% Save all data
    if runDP==0 && runP==1 && runRD==1 && runMini==0 %Standard output
        save(matfile,'-v7.3' ...
            ,'runname',...
            'COETotRD', 'RDtheory_GWh','RDPnet','Q_RD_design','Q_RD_design_LF','RDRegion_id','RDCountry_id','latOut','lonOut','ro','co','OptDH','OptDL','OptPop','OptLV','RDDepth','RDCostElements','dfhdamRD','dfldamRD','MiniFlag','OptInv', 'RDP', 'RDVolumeLake15s', 'RDSurfaceLake15s','Zoutlets', 'DisRoadOutlet', 'DisTransOutlet','DisSettlementOutlet','OptSpecCap','PopDisplacedOpt','RDlakeSurface','RDVolumeLake' ... % RD vars
            ,'n_DPinlets','aPinlet_windowMinEnd','PCOEend','PPnetend', 'PPtheoryend','Platminend','Plonminend','rPinMinend','cPinMinend','aPInletminEnd','CostElementsPMinEnd','aPinlet_windowMin','OptInvP','OptInvPMinend','PP_W','PPMinend','Pinlet_winend','dfQPminEnd','dfPLminEnd','HeadPminend','HeadraceLPminend','QDesignPinletMinend','QDesignLFPinletMinend','QDesignMeanPinletMinend','ZPinletMinend','nPipePMinend','OptSpecCapPMinend','rPin','cPin','OpthfPminend', 'OptDminend'... % P vars
            ,'PID','COEAll','PnetAll','SysID','ro_arc', 'co_arc','ros','cos','COEAlls','PnetAlls','SysIDs','ross','coss','DeselectedSites_unSort','DeselectedSites','lats','lons','latss','lonss','LakeSurfacesAlls','accAlls' ... % General vars
            ,'isfile_outlets','nd','do','do_trib_km','do_main_km','n_ielevations','landval_set','cost_lim','cost_constr','protarea_constr','navi_constr','rivermouth_constr','outletdeselect',...
            'varyRPDP_Qdesign','enablesmallPdesignQH','enablesmallPcost','varyRiverLevel_HPtype','betaC','betaP','betaCP','runDP','runRD','runP','depth_cutoff','ASYMRD','ASYMDP','dowin', 'inlet_main_sradius_km', 'inlet_trib_sradius_km', 'deselect_mainstream','inlet_win','bighydro_cutoff', 'MiniDamMax', 'deselect_mainstream', 'MiniHydro_deselect', 'MiniHydro_select', 'MiniHydro_special', 'ExistingDams_constraint', 'mangrove_constr', 'runMini'); %% Settings
        save(matfileCOEPOT,'-v7.3' ...
            ,'PID','COEAlls','PnetAlls','SysIDs','RDDepth','RDCountry_id','RDRegion_id','ro_arc', 'co_arc','ross','coss','latss','lonss'); %Just COE and Pnet
        
    elseif runDP==1 && runP==1 && runRD==1 && runMini==0
        save(matfile,'-v7.3' ...
            ,'COEend','DPPnetend','lat','lon','DPRegion_id','DPCountry_id','rin','cin','rinMinend','cinMinend','lonminend','latminend','aInletminEnd','ainlet_windowMinEnd','ahdammin','aldammin','CostElementsMinend','ainlet_windowMin','rinMin','cinMin','DPDepth' ... % DP vars
            ,'COETotRD','RDPnet','Q_RD_design','Q_RD_design_LF','RDRegion_id','RDCountry_id','latOut','lonOut','ro','co','OptDH','OptDL','RDDepth','RDCostElements','dfhdamRD','dfldamRD','MiniFlag','OptInv', 'RDP', 'RDVolumeLake15s', 'RDSurfaceLake15s','Zoutlets','DisOutlet','OptSpecCap' ... % RD vars
            ,'aPinlet_windowMinEnd','PCOEend','PPnetend','Platminend','Plonminend','rPinMinend','cPinMinend','aPInletminEnd','CostElementsPMinEnd','aPinlet_windowMin','OptInvP','OptInvPMinend','PP','PPMinend','Pinlet_winend','dfQPminEnd','dfPLminEnd','HeadPminend','QDesignPinletMinend','QDesignLFPinletMinend','QDesignMeanPinletMinend','ZPinletMinend','nPipePMinend' ... % P vars
            ,'COEAll','PnetAll','SysID','ros','cos','COEAlls','PnetAlls','SysIDs','ross','coss','DeselectedSites_unSort','DeselectedSites','lats','lons','latss','lonss' ... % General vars
            ,'isfile_outlets','nd','do','ni','discost_set','landval_set','cost_lim','cost_constr','protarea_constr','navi_constr','rivermouth_constr','outletdeselect','betaC','betaP','betaCP','runDP','runRD','runP','depth_cutoff','ASYMRD','ASYMDP','dowin', 'sradius', 'deselect_mainstream','inlet_win','inlet_window','bighydro_cutoff', 'MiniDamMax', 'deselect_mainstream', 'MiniHydro_deselect', 'MiniHydro_select', 'MiniHydro_special', 'ExistingDams_constraint', 'mangrove_constr', 'runMini'); %% Settings
        
    elseif runDP==1 && runP==0 && runRD==1 && runMini==0
        save(matfile,'-v7.3' ...
            ,'COEend','DPPnetend','lat','lon','DPRegion_id','DPCountry_id','rin','cin','rinMinend','cinMinend','lonminend','latminend','aInletminEnd','ainlet_windowMinEnd','ahdammin','aldammin','CostElementsMinend','ainlet_windowMin','rinMin','cinMin','DPDepth' ... % DP vars
            ,'COETotRD','RDPnet','Q_RD_design','Q_RD_design_LF','RDRegion_id','RDCountry_id','latOut','lonOut','ro','co','OptDH','OptDL','RDDepth','RDCostElements','dfhdamRD','dfldamRD','MiniFlag','OptInv', 'RDP', 'RDVolumeLake15s', 'RDSurfaceLake15s','Zoutlets','DisOutlet','OptSpecCap' ... % RD vars
            ,'COEAll','PnetAll','SysID','ros','cos','COEAlls','PnetAlls','SysIDs','ross','coss','DeselectedSites_unSort','DeselectedSites','lats','lons','latss','lonss','BIDs','CIDs' ... % General vars
            ,'isfile_outlets','nd','do','ni','discost_set','landval_set','cost_lim','cost_constr','protarea_constr','navi_constr','rivermouth_constr','outletdeselect','betaC','betaP','betaCP','runDP','runRD','runP','depth_cutoff','ASYMRD','ASYMDP','dowin', 'sradius', 'deselect_mainstream','inlet_win','inlet_window','bighydro_cutoff', 'MiniDamMax', 'deselect_mainstream', 'MiniHydro_deselect', 'MiniHydro_select', 'MiniHydro_special', 'ExistingDams_constraint', 'mangrove_constr', 'runMini'); %% Settings
    elseif runDP==0 && runP==0 && runRD==1 && runMini==0
        save(matfile,'-v7.3' ...
            ,'COETotRD','RDPnet','Q_RD_design','Q_RD_design_LF','RDRegion_id','RDCountry_id','latOut','lonOut','ro','co','OptDH','OptDL','RDDepth','RDCostElements','dfhdamRD','dfldamRD','MiniFlag','OptInv', 'RDP', 'RDVolumeLake15s', 'RDSurfaceLake15s','Zoutlets','DisOutlet','OptSpecCap' ... % RD vars
            ,'COEAll','PnetAll','SysID','ros','cos','COEAlls','PnetAlls','SysIDs','ross','coss','DeselectedSites_unSort','DeselectedSites','lats','lons','latss','lonss' ... % General vars
            ,'isfile_outlets','nd','do','ni','discost_set','landval_set','cost_lim','cost_constr','protarea_constr','navi_constr','rivermouth_constr','outletdeselect','betaC','betaP','betaCP','runDP','runRD','runP','depth_cutoff','ASYMRD','ASYMDP','dowin', 'sradius', 'deselect_mainstream','inlet_win','inlet_window','bighydro_cutoff', 'MiniDamMax', 'deselect_mainstream', 'MiniHydro_deselect', 'MiniHydro_select', 'MiniHydro_special', 'ExistingDams_constraint', 'mangrove_constr', 'runMini'); %% Settings
    elseif runDP==0 && runP==1 && runRD==0 && runMini==0
        save(matfile,'-v7.3' ...
            ,'latOut','lonOut','ro','co','aPinlet_windowMinEnd','PCOEend','PPnetend','Platminend','Plonminend','rPinMinend','cPinMinend','aPInletminEnd','CostElementsPMinEnd','aPinlet_windowMin','OptInvP','PP','Pinlet_winend','dfQPminEnd','dfPLminEnd','HeadPminend','QDesignPinletMinend','QDesignLFPinletMinend' ... % P vars
            ,'nd','do','ni','discost_set','landval_set','cost_lim','cost_constr','protarea_constr','navi_constr','rivermouth_constr','outletdeselect','betaC','betaP','betaCP','runDP','runRD','runP','depth_cutoff','ASYMRD','ASYMDP','dowin', 'sradius', 'deselect_mainstream','inlet_win','bighydro_cutoff', 'MiniDamMax', 'deselect_mainstream', 'MiniHydro_deselect', 'MiniHydro_select', 'MiniHydro_special', 'ExistingDams_constraint', 'mangrove_constr', 'runMini'); %% Settings
    elseif runDP==0 && runP==0 && runRD==1 && runMini==1
        save(matfile,'-v7.3' ...
            ,'COETotRD','RDPnet','Q_RD_design','Q_RD_design_LF','RDRegion_id','RDCountry_id','latOut','lonOut','ro','co','OptDH','OptDL','RDDepth','RDCostElements','dfhdamRD','dfldamRD','MiniFlag','OptInv', 'RDP', 'RDVolumeLake15s', 'RDSurfaceLake15s','OptSpecCap' ... % RD vars
            ,'isfile_outlets','nd','do','ni','discost_set','landval_set','cost_lim','cost_constr','protarea_constr','navi_constr','rivermouth_constr','outletdeselect','betaC','betaP','betaCP','runDP','runRD','runP','depth_cutoff','ASYMRD','ASYMDP','dowin', 'sradius', 'deselect_mainstream','inlet_win','bighydro_cutoff', 'MiniDamMax', 'deselect_mainstream', 'MiniHydro_deselect', 'MiniHydro_select', 'MiniHydro_special', 'ExistingDams_constraint', 'mangrove_constr', 'runMini'); %% Settings
    end
    
    
    
    disp('Data saved')
    
end
toc
fprintf('HYDRUS02 END RUN for %s: \n %s\n',runname,datetime('now'));
disp("*******************************EOF**********************************")

if savedata;  diary off; end
end