% Prepare future pop and energy data
% Pop and energy security requirements are from Wouters maps
% Energy consumption is IGCEP data compiled for my research plan
clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

popWouter=readtable(fullfile(rootf,"data\data_prep\Wouter_PopFuture2023" ,"Hist_FutPop.csv"));
energySecReq_TWh=readtable(fullfile(rootf,"data\data_prep\EnergyDemandSupply" ,"EnergyReqs_Subbasin.xlsx"),'Sheet','SubbasinEnergyinTWh');
energyreq_MWhpercapita=0.6; %MWh/capita --not convert to TWh so out pop is in millions

% In csv files, Swat is in row 1 and Kabul in row 2. for robustness, check
% basinname order in excel matches that in matlab
popWouter=popWouter([2 1 3:10],[2 1 3:24]);
energySecReq_TWh=energySecReq_TWh([2 1 3:10],:);

yrsWouter=[2015	2020:10:2080];
woutersspnames=["Prosperous", "BAU", "Downhill"]; %ssp1,2,3

for ii=1:3
    % Separate pop by ssp
    futpopSubbasin(:,:,ii)=popWouter{:, contains(popWouter.Properties.VariableNames, ["baseline",woutersspnames{ii}])};
    futpopIndus(:,ii)= futpopSubbasin(9,:,ii)';
    futpopUIB(:,ii)=sum(futpopSubbasin(1:8,:,ii))';
    % Separate energy by ssp
    energySecReq_subbasin_TWh(:,:,ii) =energySecReq_TWh{1:10, contains(energySecReq_TWh.Properties.VariableNames, ["baseline",woutersspnames{ii}])};
    energyReqUIB_TWh(:,ii)=energySecReq_subbasin_TWh(9,:,ii)';
    energyReqIndus_TWh(:,ii)=energySecReq_subbasin_TWh(10,:,ii)';
end

% Load Pak consumption vs projection
energyConsumPak_TWh=readtable(fullfile(rootf,"data\data_prep\EnergyDemandSupply" ,"EnergyReqs_Subbasin.xlsx"),'Sheet','PakDemandForecasts','DataRange','A4:E40','VariableNamesRange','A3:E3');
disp("Loaded pop and energy use files")

save(fullfile(rootf,'data','UI','data','futDemand_pop_energy.mat'),...
    "woutersspnames","yrsWouter",'futpopSubbasin', 'futpopUIB','futpopIndus',...
    'energySecReq_subbasin_TWh','energyReqUIB_TWh' , 'energyReqIndus_TWh', 'energyConsumPak_TWh')
