%% Compare sub-basin wise agricultural production based on yield and cropfrac data
% Loads three different yield data (Wouter_cfrac*yield, Wouter_yield and David_yield)

clear all
root_data = fullfile(pwd,'data');
cellsz_m=500;
% get paths to loaded data folders
root_matdata = fullfile(root_data,'UI' ,'data');
CropPrice_fname = fullfile(root_data,'data_prep','CommodityPrices','WorldBank_Commodity.xlsx'); %with additional price data from FAOSTAT

disp('Read crop yield and price data')
celly=load(fullfile(root_matdata,'CropsPot_cell_g_m2_5minLPjML.mat'));
harvest=load(fullfile(root_matdata,'CropsPot_harvest_g_m2_5minLPjML.mat'));
cropnames= readtable(fullfile(root_data,'data_prep','LPJML_Wouter', 'croptypes_SD.txt'));
davidyield=load(fullfile(root_matdata,'CropsPot.mat')); %only 1 yr data

%[WBprices,WBheaders,~] = xlsread(CropPrice_fname,'A'); % Prices of crops from World bank

%% Load subbasin
fname=fullfile(root_matdata,'UI500m_ArcGIS.mat');
load(fname,'catchments','basinlabels', 'channel','outside', 'basinmask')

%% Find 30 yr means
disp('Eval LT mean')
for i=1:32
    celly.Crops_mean(:,:,i)=mean(celly.Crops_data{i},3);
    harvest.Crops_mean(:,:,i)=mean(harvest.Crops_data{i},3);  
end

for i=1:18
    davidyield.Crops_mean(:,:,i)=davidyield.Crops_data{i};   %only one layer but aya
end
% % Plot means - heavy!!
figure(101)
sgtitle('Rainfed')
figure(102)
for i=1:16
    figure(101),subplot(4,8,i)
    imagescnan(Crops_mean(:,:,i))
    title(Crops_name{i})
end


%% Get subbasin-wise sums of production - cell based
idata=celly.Crops_mean;
ncrop=size(idata,3);
nsubbas=length(basinlabels.basinIDs)-4; % skip last 4 basins
subTot=zeros(ncrop,nsubbas);

for cropt=1:ncrop
    cdata=idata(:,:,cropt);
    for subbast=1:nsubbas
        subTot(cropt, subbast) = nansum(cdata(catchments == basinlabels.basinIDs(subbast)),'all')*cellsz_m^2;
    end
end
%
%subplot = @(m,n,p) subtightplot (m, n, p, [0.003 0.02]);
cmap=linspecer(nsubbas,'qualitative');
figure
subplot(3,1,1)
bar(subTot,'stacked') % for colorful bars
ylabel('Total gm/yr')
xticks(1:ncrop)
xticklabels(cropnames.CropNames)
xtickangle(45)
legend(basinlabels.basinnames)
applymyplotformat('Cfrac*harvest*cell area',cmap)

%% Get subbasin-wise sums of production - harvest based
idata=harvest.Crops_mean;
ncrop=size(idata,3);
nsubbas=length(basinlabels.basinIDs)-4; % skip last 4 basins
subTot=zeros(ncrop,nsubbas);

for cropt=1:ncrop
    cdata=idata(:,:,cropt);
    for subbast=1:nsubbas
        subTot(cropt, subbast) = nansum(cdata(catchments == basinlabels.basinIDs(subbast)),'all')*cellsz_m^2;
    end
end
%
%subplot = @(m,n,p) subtightplot (m, n, p, [0.003 0.02]);
cmap=linspecer(nsubbas,'qualitative');
%figure
subplot(3,1,2)
bar(subTot,'stacked') % for colorful bars
ylabel('Total gm/yr')
xticks(1:ncrop)
xticklabels(cropnames.CropNames)
xtickangle(45)
legend(basinlabels.basinnames)
applymyplotformat('harvest*cell area',cmap)

%% Get subbasin-wise sums of production - David
david_crops={
'Rainfed cereals'
'Rainfed rice'
'Rainfed maize'
'Rainfed tropical cereals'
'Rainfed pulses'
'Rainfed roots and tubers'
'Rainfed oil crops'
'Biofuel sugarcane'
'Biofuel maize'
'Biofuel woody'
'Biofuel non-woody'
'Irrigated temperate cereals'
'Irrigated rice'
'Irrigated maize'
'Irrigated tropical cereals'
'Irrigated pulses'
'Irrigated roots and tubers'
'Irrigated oil crops'};

idata=davidyield.Crops_mean;
ncrop=18;
nsubbas=length(basinlabels.basinIDs)-4; % skip last 4 basins
subTot=zeros(ncrop,nsubbas);

for cropt=1:ncrop
    cdata=idata(:,:,cropt);
    for subbast=1:nsubbas
        subTot(cropt, subbast) = nansum(cdata(catchments == basinlabels.basinIDs(subbast)),'all')*cellsz_m^2;
    end
end
%
%subplot = @(m,n,p) subtightplot (m, n, p, [0.003 0.02]);
cmap=linspecer(nsubbas,'qualitative');
%figure
subplot(3,1,3)
bar(subTot,'stacked') % for colorful bars
ylabel('Total gm/yr')
xticks(1:ncrop)
xticklabels(david_crops)
xtickangle(45)
legend(basinlabels.basinnames)
applymyplotformat('David_harvest*cell area',cmap)