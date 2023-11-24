%% Evaluate agricultural yield based on Lpjml outputs and get ag land value
% decided to use sum of all crops rather than using just one crop, that is
% done in another code. here trying David's version based on max crop
% Based on Yield_maker.m by David

clear all
root_data = fullfile(pwd,'data');
cellsz_m=500;
% get paths to loaded data folders
root_matdata = fullfile(root_data,'UI' ,'data');
CropPrice_fname = fullfile(root_data,'data_prep','CommodityPrices','WorldBank_Commodity.xlsx'); %with additional price data from FAOSTAT

selYears=20:30; % data is from 1981-2010, select 2001-2010 because costs are higher at the end of the period
%% Redo David's method w Wouter data
harvest=load(fullfile(root_matdata,'CropsPot_harvest_g_m2_5minLPjML.mat')); %yield only
davidyield=load(fullfile(root_matdata,'CropsPot_David.mat')); %only 1 yr data
load(fullfile(root_matdata,'CropsPot_Lpjml.mat'),'pCrop_USDPerTonne','cropnames_sel');
disp('Read crop yield and price data')
tmp=strcat('Irri: ',cropnames_sel(:) );
Crops_name={cropnames_sel{:} ,' ' ,' ' ,' ' , tmp{:},' ' ,' ' ,' ' };
pCrop_USDPerTonne=[pCrop_USDPerTonne 0 0 0 pCrop_USDPerTonne 0 0 0];

%% Crops value map as product of rainfed crop production and price
for i=1:32   
    tmp=harvest.Crops_data{i};
    Crops_data(:,:,i)= mean(tmp(:,:,selYears),3);
    Crops_value(:,:,i) = Crops_data(:,:,i) * pCrop_USDPerTonne(i); % t/km^2 * $/t = $/km^2
end
disp('Evaluated value per crop')


%% Mask basin
load(fullfile(root_matdata,'UI500m_ArcGIS.mat'), 'basinmask')
Crops_value= maskBasin(Crops_value,basinmask);

%% Select max of 13*2 rainfed crop type value per cell
selcrops=[1:13 17:29]; % Skip the  last 3 categories
Crops_value_sel=Crops_value(:,:,selcrops);
[Crops_max,maxID]=max(Crops_value_sel,[],3);

%Count of maxID
dominantCrop=countUniques(maxID);
t=table(dominantCrop(:,2),'RowNames',Crops_name(selcrops(dominantCrop(:,1))));
disp('Build max value map')

figure
bar(dominantCrop(:,2))
xticklabels(Crops_name(selcrops(dominantCrop(:,1))))
xtickangle(45)
title("Maximum value crop in majority cells")

%% Plot land value map and spatial sum for each croptype
figure
subplot(1,2,1)
imagescnan(maskBasin(Crops_max,basinmask))
title('Maximum land value map based on potential rainfed+irrigated food production (USD/km^2)','Fontsize',12)
set(gca, 'ColorScale', 'log')
mycbar("Agricultural value (USD/km^2)")
subplot(1,2,2)
bar(squeeze(nansum(Crops_value_sel,[1 2])))
set(gca,'xtick',1:26,'xticklabel',Crops_name(selcrops),'XTickLabelRotation',45)
ylabel("Agricultural value (USD/km^2)")
grid on

%% Save calculated landvalue data
LandValue_USD_km2=Crops_max;
disp('Saving')
matfile = fullfile(root_matdata, sprintf('LandValue_2001-2010_maxLpjml.mat'));
save(matfile,'LandValue_USD_km2','-v7.3');
disp("***************************EOF******************************")
