% Prepare land acquisition cost matrix. Reclassify ESA Land Use (LU) map 
% to match Dasu land prices class. First reclassify
% 37 ESA LU classes to 5 generic classes (Urban, Water, Cultivated,
% Uncultivated, Barren). Then use HYDE population map to explore potential
% population density thresholds to divide are into rural and urban. Using
% identified threshold, split cultivated, uncultivated and barren to
% rural-urban. Finally, used settlement map to set settlement cells w a
% buffer as residential area. Classify residential rural-urban same as
% before. Final classification is used to generate cost map that is saved
% in LandAcq.mat file

% Refer to LULC classes excel sheet for data
% Created 12 Nov 2020

clear all
close all

% Load  basin bounds
load('G:\SurfDrive\GitConnect\data\UI\data\Basin\Basin_551.mat','outside')

%% Load ESA LU map
load('G:\SurfDrive\GitConnect\data\UI\data\Landuse.mat')
ESA_LC=maskBasin(data,~outside);
figure;imagescnan(ESA_LC)
axis image

%% Reclassify 37 ESA LU map to 5 Dasu land classes
OldLU=[nan 0	10	11	12	20	30	40	50	60	61	62	70	71	72	80	81	82	90	100	110	120	121	122	130	140	150	152	153	160	170	180	190	200	201	202	210	220];
Old2NewLU=[nan 0	2	3	3	2	2	2	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	4	4	4	4	5	5	5	1	4	4	4	5	5];
NewLU=1:5;
NewLUnames= {'Urban', 'Cultivated', 'Uncultivated', 'Barren', 'Water/Snow'};

Dasu_LC=changem(ESA_LC,Old2NewLU,OldLU);
figure;imagescnan(Dasu_LC)

%check classification 
unique(ESA_LC(Dasu_LC==5))
figure;imagescnan(Dasu_LC==1)
title("Urban grids in ESA")
save('G:\SurfDrive\GitConnect\data\UI\data\Landuse_Dasu1.mat','Dasu_LC')

%% Load HYDE population and discretize into rural-urban
load('G:\SurfDrive\GitConnect\data\UI\data\Popd.mat')
popUI=maskBasin(data,~outside);
figure; imagescnan(popUI<500& popUI>0)

% Plots of population distribution
figure; plot(sort(popUI(:)))
set(gca, 'Yscale', 'log')
applymyplotformat()
ylabel("Population density (inhabitants/km2)")
%
figure
cdfplot(popUI(:))
xlabel("Population density (inhabitants/km2)")
set(gca, 'Xscale', 'log')
title('Population density across UI')


%% Load settlements
load('G:\SurfDrive\GitConnect\data\UI\data\Settlements.mat')
settlement=maskBasin(data,~outside);
[settlement_r, settlement_c]=find(maskBasin(data,~outside)>0);
settlement_i=find(maskBasin(data,~outside)>0);
s_pop=popUI(settlement_i);

%
figure
cdfplot(s_pop(:))
xlabel("Population density (inhabitants/km2)")
set(gca, 'Xscale', 'log')
title('Population density in settlements in UI')

%% Create rural urban class based on population
ruralthres=[100 200 250 500]; %upper limit in population per km2
popUI_rural500 = discretize(popUI, [0 500 Inf]); % ,'categorical',{'R', 'U'}); 
figure;
for i=1:4
    subplot(2,2,i)
    poptmp=discretize(popUI, [0 ruralthres(i) Inf]);
    imagescnan(poptmp)
    hold all
    scatter(settlement_c, settlement_r,rescale(s_pop,2,100),'ko')
    title(sprintf("Rural is <= %d", ruralthres(i)))
    
  % Compare ESA urban and pop density urban
    fprintf("No of urban cells overlapping w ESA for %d threshold: %d\n", ruralthres(i), numel(find((Dasu_LC==1) & (poptmp==2))))
    %figure;imagescnan(maskBasin(Dasu_LC==1 & pop250==2,~outside))
end
sgtitle('Blue=Rural, Yellow=Urban')
legend('settlement by pop density')


%% Reclassify Dasu LULC into rural and urban
pop250=discretize(popUI, [0 ruralthres(3) Inf]);
NewLU2=101:109;
NewLUnames2= {'Residential - Urban', 'Residential - Rural', 'Cultivated - Urban', 'Cultivated - Rural', 'Uncultivated - Urban', 'Uncultivated - Rural', 'Barren - Urban', 'Barren - Rural', 'Water/Snow'};

%initialize
Dasu_LC2=Dasu_LC*0+100;

figure
colormap(parula(10))
figc=1;
subplot(2,5,figc)
imagescnan(Dasu_LC)
title("Old Dasu")
colorbar
% figc=2;
% subplot(2,5,figc)
% %imagescnan(Dasu_LC2)
% title(NewLUnames2(figc))
% colorbar

% classify cultivated, uncultivated and barren to urban and rural
for c=2:4
    Dasu_LC2(Dasu_LC==c & pop250 == 2)=NewLU2(c*2-1); %urban
    subplot(2,5,figc+1)
    imagescnan(Dasu_LC2)
    title(NewLUnames2(c*2-1))
    colorbar
    
    Dasu_LC2(Dasu_LC==c & pop250 == 1)=NewLU2(c*2); %rural
    subplot(2,5,figc+2)
    imagescnan(Dasu_LC2)
    title(NewLUnames2(c*2))
    figc=figc+2;
    colorbar
end

%% Checks on new LU data
disp('Checks')
unique(Dasu_LC2(~isnan(Dasu_LC2)))'
all(all((Dasu_LC2==104|Dasu_LC2==106|Dasu_LC2==108) == (pop250==1 & Dasu_LC~=5 & Dasu_LC~=1)))  % all rural should be same
all(all((Dasu_LC2==101|Dasu_LC2==103|Dasu_LC2==105|Dasu_LC2==107) == (pop250==2 & Dasu_LC~=5)))  % all urban should be same
% find indices that dont match in urban classification
chkidx=(Dasu_LC2==101|Dasu_LC2==103|Dasu_LC2==105|Dasu_LC2==107)-(pop250==2 & Dasu_LC~=5 );
unique(Dasu_LC(chkidx==-1)) % shows that all are urban
unique(pop250(chkidx==-1)) % shows that all are rural...so these are cells that are rural in pop map and urban in ESA!! this is fine for now

%% Adding residential: for settlements
buffer_settlement=2; %2*500m=1km which means each settlement is assumed to be spread over 9km2
s_ruralurban=settlement*0;
s_ruralurban(settlement_i)=pop250(settlement_i);
% identify settlements that are in urban vs rural areas and add buffer
% around them
ruralrescells=logical(createBuffer(s_ruralurban==1,buffer_settlement));
urbanrescells=logical(createBuffer(s_ruralurban==2,buffer_settlement));

% check what land use these are currently in -some rural cells are in urban
disp('Count of LUclasses in settlement cells and buffers')
countUniques(Dasu_LC2(s_ruralurban==1)) 
countUniques(Dasu_LC2(s_ruralurban==2))
countUniques(Dasu_LC2(ruralrescells))
countUniques(Dasu_LC2(urbanrescells))

% Reclassify settlement cells into urban or rural residential
Dasu_LC2(ruralrescells)=NewLU2(2);
subplot(2,5,figc+1)
imagescnan(Dasu_LC2)
title(NewLUnames2(2))
colorbar
%% Adding urban residential: from settlement and ESA urban
Dasu_LC2(urbanrescells)=NewLU2(1);
Dasu_LC2(Dasu_LC==1)=NewLU2(1);
subplot(2,5,figc+2)
imagescnan(Dasu_LC2)
title(NewLUnames2(1))
colorbar

%% Convert ESA water to water
Dasu_LC2(Dasu_LC==5)=NewLU2(end);
subplot(2,5,figc+3)
imagescnan(Dasu_LC2)
title(NewLUnames2(end))
colorbar


%% Evaluate land acquisition cost using land rates for Dasu Hydropower Project
NewLU2=101:109;
NewLU2_USDperkm2 = [23203800.44,	16875491.23,	23203800.44,	16875491.23,	8437745.61,	4858321.33,	5265153.26,	2109436.40, 0];

LandClass=Dasu_LC2;
LandClass(outside)=0;
disp('Final count of Land Use classes:')
countUniques(LandClass)

LandAcqVal_USDperkm2=zeros(size(Dasu_LC2));

for i=1:length(NewLU2)
    LandAcqVal_USDperkm2(LandClass==NewLU2(i))= NewLU2_USDperkm2(i);
end

%% Final figures
figure
imagescnan(LandAcqVal_USDperkm2)
title("Land acquisition cost in USD per km^2")
colorbar

figure
imagescnan(ESA_LC)
title("ESA LU classes")
colorbar

figure
imagescnan(Dasu_LC)
title("ESA LU classes converted to 5 classes")
colormap(parula(5))
ch=colorbar;
ch.Ticks=NewLU;
ch.TickLabels=NewLUnames;


figure
imagescnan(Dasu_LC2)
title("Final land classes considering settlements and rural-urban population")
colormap(parula(9))
ch=colorbar;
ch.Ticks=NewLU2;
ch.TickLabels=NewLUnames2;

%% Save outputs
save('G:\SurfDrive\GitConnect\data\UI\data\LandAcq.mat', 'LandClass', 'LandAcqVal')

load('G:\SurfDrive\GitConnect\data\UI\data\Basin\Basin_551.mat','Rw','proj')
geotiffwrite('G:\SurfDrive\GitConnect\data\UI\data\LandAcq_LUclasses.tif',LandClass,Rw,'GeoKeyDirectoryTag',proj.GeoTIFFTags.GeoKeyDirectoryTag)
geotiffwrite('G:\SurfDrive\GitConnect\data\UI\data\LandAcq_USDperkm2.tif',LandAcqVal_USDperkm2,Rw,'GeoKeyDirectoryTag',proj.GeoTIFFTags.GeoKeyDirectoryTag)
