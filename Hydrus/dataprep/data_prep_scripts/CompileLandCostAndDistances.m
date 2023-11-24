%% Compile the three land cost classes and 3 distance classes
clear all
close all

% Load land cost and distance mats
load(fullfile(pwd,'data\UI\data\LandAcq.mat'), 'LandAcqVal')
load(fullfile(pwd,'data\UI\data\LandValue_2001-2010_sumLpjml.mat'), 'LandValue_USD_km2')
load(fullfile(pwd,'data\UI\data\TreeCost.mat'))
Costdatanames={'Tree Loss', 'Agricultural Loss', 'Land Acquisition'};

% Convert distances to km
cellsz=0.500; % in km
load(fullfile(pwd,'data\UI\data\Dis2Road.mat'), 'data')
Dis2Road=data*cellsz;
load(fullfile(pwd,'data\UI\data\Dis2Transmission.mat'), 'data')
Dis2Transmission=data*cellsz;
load(fullfile(pwd,'data\UI\data\Dis2Settlement.mat'), 'data')
Dis2Settlement=data*cellsz;

%% nan out area outside the basin
load(fullfile(pwd,'data\UI\data\Basin\Basin_551.mat'),'outside')
%Load new basin
%load(fullfile(pwd,'data\UI\data\UI500m_ArcGIS.mat', 'outside')

TreeCost(outside)=nan;
LandAcqVal(outside)=nan;
LandValue_USD_km2(outside)=nan;
Dis2Road(outside)=nan;
Dis2Transmission(outside)=nan;
Dis2Settlement(outside)=nan;

%% Plot 3 distances as subfigure
subplot = @(m,n,p) subtightplot (m, n, p, [0.003 0.02]);

figure;
subplot(2,3,1)
plotMyRaster(Dis2Road, "Nearest road (km)")

subplot(2,3,2)
plotMyRaster(Dis2Transmission, "Nearest transmission line (km)")

subplot(2,3,3)
plotMyRaster(Dis2Settlement, "Nearest settlement (km)")
%sgtitle("Distance constraints on technical and economic potential")

%% Plot 3 land value as subfigure
%figure
subplot(2,3,4)
plotMyRaster(LandAcqVal,  "Acquisition value (USD/km^2)",[],1)

subplot(2,3,5)
plotMyRaster(TreeCost,  "Tree value (USD/km^2)",[],1)

subplot(2,3,6)
plotMyRaster(LandValue_USD_km2, "Agricultural value (USD/km^2)",[],1)

%sgtitle("Land use constraints on technical and economic potential")
%% For saving as PDF
set(gcf, 'Color', 'w');
export_fig(fullfile(pwd,'data\UI\CompiledDistance+CostMaps5.pdf'))
export_fig('G:\OneDrive - ICIMOD\PhD_Paper1\Paper1_Fig\R0_Figs\8_LandValueMaps_R2.jpg')
export_fig('G:\OneDrive - ICIMOD\PhD_Paper1\Paper1_Fig\R0_Figs\8_LandValueMaps_R2.pdf')

orient('landscape')
print(fullfile(pwd,'data\UI\CompiledDistance+CostMaps3.pdf','-dpdf','-bestfit'))

%
export_fig('G:\OneDrive - ICIMOD\PhD_Paper1\Paper1_Fig\R0_Figs\CompiledDistance.pdf')
export_fig('G:\OneDrive - ICIMOD\PhD_Paper1\Paper1_Fig\R0_Figs\CostMaps.pdf')

%% Compare 3 land associated costs
figure
plot(TreeCost(:),'.','Color',[0.2314    0.4392    0.0784])
hold all
plot(LandValue_USD_km2(:),'o','Color',[0.8902    0.2118    0.2118])
plot(LandAcqVal(:),'x','Color',[0.9294    0.6941    0.1255])
set(gca, 'Yscale', 'log')
grid on
xlabel('Cell ID')
ylabel('Cost in USD per km^2')
legend(Costdatanames)

%% Compare 3 land associated costs, sorted
figure
plot(sort(TreeCost(:)),'.','Color',[0.2314    0.4392    0.0784])
hold all
plot(sort(LandValue_USD_km2(:)),'o','Color',[0.8902    0.2118    0.2118])
plot(sort(LandAcqVal(:)),'+','Color',[0.9294    0.6941    0.1255])
set(gca, 'Yscale', 'log')
grid on
xlabel('Sorted cells')
ylabel('Cost in USD per km^2')
legend(Costdatanames)

%% Cross plots
figure
plot(LandAcqVal(:),LandValue_USD_km2(:),'o')
hold all
%figure
plot(LandAcqVal(:),TreeCost(:),'+')
grid on