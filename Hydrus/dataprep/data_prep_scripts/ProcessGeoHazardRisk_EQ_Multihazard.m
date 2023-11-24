% Compare EQ data w 3 levels or 5 levels of hazard. Prepped data is saved
% and added to basin compilation
% Check multi hazard risk levels to figure out how to represent single and
% multi-hazard classes. The multi-hazard score is assessed in the runHydrus
% code itself. so here it is just trial, saved data is bit used later.

clear
close all

buffer_thrust_km = 1;       %Buffer for thrust in km
buffer_fault_km = 1;        %Buffer for fault in km
buffer_glof_km = 0.5;       %Buffer for glof paths in km

geohazard_select    =[1 1 1 1]; %Apply EQ-PGA, EQ-Thrust, LS, GLOF hazard risks
geohazard_scen      =1;
geohazard_level_rates =[
    0	0  % Added cost for each hazard level 1=Very Low, 5=Very High
    1	2
    2	4
    3	6
    4	8
    5	10];
old3levels=0;
new3levels=0;

hazard_lnames={'Very Low', 'Low', 'Medium', 'High','Very High'};
outfname = fullfile(pwd, 'data', 'UI', 'data', 'geohazards_compiled.mat');

%% oldGeoHazard 3 classes
if old3levels
    basinfname = fullfile(pwd, 'data', 'ASIA', 'Basin_UIB', sprintf('PantpeBasin_%d.mat', 3));
    load(basinfname,'landslide_level','seismic_level','seismicthrust', 'glofpath','cellsz_m','outside')
    
    geohazard_highrisk = (seismic_level==3) + createBuffer(seismicthrust==1, buffer_thrust_km*1e3/cellsz_m)...
        + createBuffer(seismicthrust>=2, buffer_fault_km*1e3/cellsz_m)...
        + (landslide_level==3) +createBuffer(glofpath,buffer_glof_km*1e3/cellsz_m)*10;
    
    geohazard_medrisk =(seismic_level==2)+(landslide_level==2);
    
    
    figure;subplot(1,2,1);imagescnan(geohazard_highrisk);
    title(sprintf('High geo hazard risk: %0.2f%% of basin', sum(geohazard_highrisk>0,'all')/sum(~outside,'all')*100))
    subplot(1,2,2);imagescnan(geohazard_medrisk);
    title(   sprintf('Med geo hazard risk: %0.2f%% of basin', sum(geohazard_medrisk>0,'all')/sum(~outside,'all')*100))
    sgtitle("Old LS 3 hazard level data")
end

%% Load new LS and use as 3 level hazard
basinfname = fullfile(pwd, 'data', 'ASIA', 'Basin_UIB', sprintf('PantpeBasin_%d.mat', 4));
load(basinfname,'landslide_level','seismic_level','seismicthrust', 'glofpath','outside','cellsz_m')
nbasin_cells=sum(~outside,'all');
landslide_5level=landslide_level;
seismic_3level=seismic_level;

if new3levels
    landslide_3level=changem(landslide_5level,[0 1 1 2 2 3],0:5); % Reclass LS data
    
    classname={'Low','Med','High'};
    oldEQclass=0;
    if oldEQclass
        
        seismic_3level=reclassify(seismicpga_data,[1.5, 4]*g, 1:3, classname);
        %all(seismic_level2(~isnan(seismic_level2(:)))==seismic_level(~isnan(seismic_level(:))),'all')
    end
    
    geohazard_highrisk = (seismic_3level==3) + createBuffer(seismicthrust==1, buffer_thrust_km*1e3/cellsz_m)...
        + createBuffer(seismicthrust>=2, buffer_fault_km*1e3/cellsz_m)...
        + (landslide_3level==3) +createBuffer(glofpath,buffer_glof_km*1e3/cellsz_m)*10;
    
    geohazard_medrisk =(seismic_3level==2)+(landslide_3level==2);
    
    figure;subplot(1,2,1);imagescnan(geohazard_highrisk);
    title(sprintf('High geo hazard risk: %0.2f%% of basin', sum(geohazard_highrisk>0,'all')/sum(~outside,'all')*100))
    subplot(1,2,2);imagescnan(geohazard_medrisk);
    title(   sprintf('Med geo hazard risk: %0.2f%% of basin', sum(geohazard_medrisk>0,'all')/sum(~outside,'all')*100))
    sgtitle('New LS - 3 hazard level')
end

%% Get 5 level EQ data
seismicpga_data = load(fullfile(pwd, 'data', 'UI\data\GSHAP.mat'));
seismicpga_data = seismicpga_data.data;

% Redo class EQ risk
g=9.81; % accelaration due to gravity
seismic_5level=reclassify(seismicpga_data,[.16 .22 .32 .408]*g, 1:5);
seismic_5level(outside) = 0;
%seismic_3level=changem(seismic_5level,[0 1 1 2 2 3],0:5); % Reclass EQ data

figure;
subplot(1,2,1)
imagescnan(seismic_3level)
subplot(1,2,2)
imagescnan(seismic_5level)
sgtitle('EQ 3 level vs 5 level')

%% Compile all hazards
geohazard_names={'Earthquake-PGA', 'Earthquake-Thrust', 'Landslide', 'GLOF'};
nhazardlevel=5;
geohazard_compile(:,:,1)=geohazard_select(1)*seismic_5level;           % seismic_5level
geohazard_compile(:,:,2)=geohazard_select(2)*createBuffer(seismicthrust>0, buffer_thrust_km*1e3/cellsz_m)*nhazardlevel;  % seismicthrust_level
geohazard_compile(:,:,3)=geohazard_select(3)*landslide_5level;          % landslide_level Reclassified
geohazard_compile(:,:,4)=geohazard_select(4)*createBuffer(glofpath,buffer_glof_km*1e3/cellsz_m)*nhazardlevel;  % glof_levelGLOF  path + buffer considered highest risk level

%% Explore multi hazard -- tried sum and occurence decided to use sum 
geohazard_sum=sum(geohazard_compile,3);
geohazard_occurence = 0*geohazard_sum;
for i=1:size(geohazard_compile,3)
    geohazard_occurence=geohazard_occurence+double(geohazard_compile(:,:,i)>0);
end

figure;
subplot(2,2,1)
imagescnan(geohazard_occurence)
colorbar
title("Number of hazards in each cell")
subplot(2,2,2)
imagescnan(geohazard_sum)
title("Sum of hazards in each cell")
colorbar
subplot(2,2,3)
histogram(geohazard_occurence(geohazard_occurence>0))
set(gca, 'YScale', 'log')
grid on
subplot(2,2,4)
histogram(geohazard_sum(geohazard_sum>0))
set(gca, 'YScale', 'log')
grid on

%% Get max hazard classes and rates
% Max hazard for Status_Quo and Risk Averse: Get hazard rate based on max hazard class -  this is kept in main code as will vary for # of hazards selected
% Composite considering max of all hazards
geohazard_max_5levels=max(geohazard_compile,[],3);
%    geohazard_max_5levels=maskBasin(max(geohazard_compile,[],3), ~outside);
geohazard_hazard_rates_max=changem(geohazard_max_5levels,geohazard_level_rates(:,2),geohazard_level_rates(:,1));

%% Get multi hazard classes and rates Composite considering sum of all hazards rescaled from 0-5
geohazard_sum=maskBasin(sum(geohazard_compile,3), ~outside);

rescalemultabs=0; % minmax in geohazard_sum is 2-19. if rescalemultabs then rescale 0-20 to 0-5 else rescale data min max from 0-5
%if rescalemultabs
    geohazard_multi_score1=rescale(geohazard_sum,0,5,'InputMin',0,'InputMax',20);
    geohazard_multi_5levels1=reclassify(geohazard_multi_score1, 1:4, 1:5); %,hazard_cnames);  % need to this as background vals are 0
    figure;
        subplot(3,2,1)
    imagescnan(geohazard_multi_score1);title("Absolute Scale 0-20 to 0-5 - SCORE")
    colorbar("Location",'southoutside')
        subplot(3,2,3)
    imagescnan(geohazard_multi_5levels1);title("Absolute Scale 0-20 to 0-5 - HAZARD LEVEL")
    colorbar("Location",'southoutside')
        subplot(3,2,5)
    histogram(geohazard_multi_5levels1(:))
    set(gca, 'YScale', 'log')
    grid on
    title("Hazard level distn")

    xticklabels(hazard_lnames)
    countUniques(geohazard_multi_5levels1(:), nbasin_cells)
%
%else
    geohazard_multi_score=rescale(geohazard_sum,0,5);%,'InputMin',2,'InputMax',19);
    %geohazard_multi_5levels=reclassify(geohazard_multi_score, 1:4,1:5);
    %%,hazard_cnames);  % both this and the one below are the same
       geohazard_multi_levels=reclassify(geohazard_multi_score, 0:4,0:5); %,hazard_cnames);

    subplot(3,2,2)
    imagescnan(geohazard_multi_score);title("Relative Scale min-max to 0-5 - SCORE")
    colorbar("Location",'southoutside')
    
    subplot(3,2,4)
    imagescnan(geohazard_multi_5levels);title("Relative Scale min-max to 0-5 - HAZARD LEVEL")
        colorbar("Location",'southoutside')

    subplot(3,2,6)
    histogram(geohazard_multi_5levels(:))
    title("Hazard level distn")
    set(gca, 'YScale', 'log')
grid on
xticklabels(hazard_lnames)
countUniques(geohazard_multi_5levels(:), nbasin_cells)
   
%end
tmp=[countUniques(geohazard_multi_5levels1(:), nbasin_cells)  countUniques(geohazard_multi_5levels(:), nbasin_cells)]

%Get rates for hazard classes
geohazard_hazard_rates_multi=changem(geohazard_multi_5levels,geohazard_level_rates(:,2),geohazard_level_rates(:,1));



%% Plot individual, max and multi hazard risks
figure
for i=1:size(geohazard_compile,3)
    %countUniques(geohazard_compile(:,:,i))
    subplot(3,4,i)
    imagescnan(geohazard_compile(:,:,i))
    title(geohazard_names(i))
    % 
end
colorbar
subplot(3,4,5:6)
imagescnan(geohazard_max_5levels)
title("Maximum hazard based risk levels")
%colorbar

subplot(3,4,7:8)
imagescnan(geohazard_multi_5levels)
title("Multi hazard based risk levels")
colorbar('Location','southoutside')

subplot(3,4,9:10)
imagescnan(geohazard_hazard_rates_max)
title("Maximum hazard based added %")
%colorbar

subplot(3,4,11:12)
imagescnan(geohazard_hazard_rates_multi)
title("Multi hazard based added %")
colorbar('Location','southoutside')

%% Cell counts for the max and multi hazard cases
mm=[countUniques(geohazard_max_5levels(geohazard_max_5levels>0), nbasin_cells) ...
    countUniques(geohazard_multi_5levels(geohazard_multi_5levels>0), nbasin_cells)];
figure
bar(mm(:,[3,6]))
xticklabels(hazard_lnames)
legend(["Maximum hazard based", "Multi hazard based"])
title ('Percentage of cells in each risk level')
grid on

%% Save new hazard classes for single and multi hazards -- the max hazard one is autogenerated in the main code
save(outfname, 'landslide_5level', 'seismic_5level', 'geohazard_max_5levels', 'geohazard_multi_5levels', ...
    'geohazard_names', 'nhazardlevel','hazard_lnames');