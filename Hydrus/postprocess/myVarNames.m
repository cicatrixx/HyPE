% Standardized parameter names, labels and colors to use across plots
addpath(genpath(fullfile(pwd,'..','GitConnect\Hydrus\devFiles\')), ...
    genpath(fullfile(pwd,'Hydrus\')))

rootf=pwd; %'G:\SurfDrive\HPmodel'; %
rootof=fullfile(rootf,'output','HistRuns_Figs_Analysis');%,'Figs_trial'

continent_in='ASIA_Historical_Runs';
nbasin=103;
costlim=0.1; %$0.01/kWh COE limit

pottypes5={'Theoretical','Technical','Financial','Sustainable','Visualized'};
pottypes5_short={'Theory','Tech','Fin','Sust','Vis'};
pottypes6=[pottypes5(:)', {'Existing'}];
pottypes3={'Technical','Financial','Sustainable'};
pottypes3_short={'Tech','Fin','Sust'};

scenarios={'Full','Remain'};
planttypes={'River Power Plant','Diversion Canal Plant'};
searchtypes={'Large', 'Medium', 'Mixed'};

tech2sust_Labels={'TECHNICAL potential'
    '+ Financial limit'
    '+ Sustainable discharge'
    '+ Sustainable areas'
    '+ Hazard risk aversion'
    'SUSTAINABLE potential'};

geohazard_scennames={'No geohazard', 'CostBased','RiskAverse', 'Multi-hazard'}; %for filenames
geohazard_scennames_cl={'No geohazard', 'Cost-based','Risk Averse', 'Multi-hazard'};% for labeling
geohazard_names={'Earthquake', 'GLOF','Landslide'};  % short ones used to gen filename
geohazard_names_short={'EQ','EQ', 'LS', 'GLOF'};  % short ones used to gen filename


% add lines and labels for different sizes of HP -- from Siddiqui
HPclass=["Mega (>1000 MW)", "Large (500-1000 MW)", "Medium (50-500 MW)", "Small (5-50 MW)", "Mini (0.15-5 MW)", "Micro (0.005-0.15 MW)", "Pico (<0.005 MW)"];
HPsz_greaterthan=[1000	500	50	10	5	0.005]/1000*365*24; %MW converted to GWh assuming year round production

%% Create color for plots
cl_RP= [216,179,101]/255; %brownish %[112 173 71]/255; %'#70AD47';
cl_DP= 127/255*[1 1 1]; %gray %'#7F7F7F';
cl_channel=  [.53 .81 1]; % very light river %[0 0.45 0.74]; %river blue

mycolors={cl_RP cl_DP};
baralpha=0.7;
myalpha=0.6; % for scatter
mygraylines=.5*[1 1 1];  % for axes borders and boxes
subbasinalpha=0.9;
% For HP classes in theory pot bar chart    colormap(cbrewer2('BrBG',length(tblsum.HPtype)))

%% subbasin colors: inferno to match wouter
%cmap8=inferno(14);
%cmap8=cmap8(4:11,:);
wouteralpha=0.8;
cmap8_inferno=[71 59 101 
    109 59 134
    145 71 136
    181 84 125
    215 103 104
    239 133 78
    250 172 53
    247 215 96]/255;

% Tick locations for subbasin colorbar
% colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'Ticklabels',(basindata.basinnames),'TickDirection','none')

%% subbasin colors: qualitative colorblind friendly combo from Tol palette
%https://personal.sron.nl/~pault/
cmap8_light=[119, 170, 221
    153,221,255
    68,187,153
    170,170,0
    238,221,136
    238,136,102
    255,170,187
    221,221,221
    ]/255;

cmap8_bright=[
    68,119,170
%136,34,85 %wine
51,34,136 %indigo
34,136,51
    204,187,68
    238,102,119
    170, 51, 119
        102,204,238
  
    187,187,187
    ]/255;

%% subbasin colors: from Wong
%from: https://www.nature.com/articles/nmeth.1618/figures/2
cmap8_wong= [
    230, 159, 0     %Orange
    86,180,233      %Sky Blue
    0, 158, 115     %Bluish green
    240, 228, 66    %Yellow
    0, 114, 178     %Blue
    213, 94, 0      %Vermillion
    0  0 0          %Black
    204,121,167     %Reddish purple
    ]/255;

cmap8_wong4R= cmap8_wong'*255;


% in AI i change the opacity to 60% for these sub-basin colors in the onion
% plot

%% Potential types onion colors
cmap6_pottype=[cbrewer2('Oranges',5); 0 0 0];

%% Energy and hazard scenarios
%cmap3_mainenergy=brighten(flipud(cbrewer2('Set2',3)),-0.4); % color blind and printer friendly purple orange green
cmap3_mainenergy=flipud(cbrewer2('Set2',3)); % color blind and printer friendly purple orange green

%% Manually assign xy for onions
%'Kabul'  'Swat'  'Indus''Jhelum''Chenab'    {'Ravi'}    {'Beas'}  {'Satluj'} 'All basins'
bas_x =[670     1050     1800    1360   1600    1400    1900    2300 2400];
bas_y =[410     380     500     600     890    1200    1100    1250 300] ;

%%
 mlaborange=[0.8500 0.3250 0.0980];
    mlabblue=[0 0.4470 0.7410];