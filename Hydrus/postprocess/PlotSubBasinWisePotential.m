%% Get sub-basin wise totals for full and remaining scenarios
% Load and compile data for main energy and hazard scenarios in pcsout
% Load visualized and theoretical data
% Get sub-basin wise total potential, per area potential and per capita potential
% Last row in sub-basin matrices have basin level totals! Cols are Theory,
% Tech, Fin, Sust, Vis and Existing potential.
% For R plotting, RP and DP type totals is only done for subbasins and not
% for existing case
clc
clearvars
close all
addpath(genpath('G:\SurfDrive\GitConnect\Hydrus\devFiles\'), ...
    genpath('G:\SurfDrive\HPmodel\Hydrus\'))
run('myVarNames.m')

ofmatname=fullfile(rootof,'Figs_trial','MainScenarios.mat');
ofxls=fullfile(rootof,'Figs_trial','SubBasinData4WouterBar.xlsx');

idx=1;
plotallfigs=0;
write2mat=1;
save4saurav=0;
%subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.06],[.2 .1],[.12 .03]); % ;
%subplot = @(m,n,p) subtightplot(m,n,p,gap [vgap hgap],marg_h [low up],marg_w [left right],varargin)
%% Load catchment and population density data
datain=load(fullfile(rootf,sprintf('\\data\\ASIA\\Basin_UIB\\PantpeBasin_%d.mat',nbasin)),'Regions', 'Pop','Countries');
cellsz_m=500;

load(fullfile(rootf,'\data\UI\data\UI500m_ArcGIS.mat'), 'basinlabels')
basindata=basinlabels(1:8,1:3);
label_subbasin_basin =[basindata.basinnames{:}, "All basin"];

% Get subbasin-wise area and population
for st=1:height(basindata)
    bidx= datain.Regions == basindata.basinIDs(st);
    basindata.Area_m2(st,1)=sum(bidx,"all")*cellsz_m^2;
    basindata.Pop_ppl(st,1) = sum(datain.Pop(bidx))*(cellsz_m/1e3)^2 ; %inhabitants /km2 * cell area in km ; %pop=population density in inhabitants/km2
end

catchments=datain.Regions;
catchments_cl=maskBasin(catchments,catchments>0&catchments<109);

% Basin area and pop
bidx= datain.Regions >0;
basinArea_m2=nansum(catchments_cl>0,"all")*cellsz_m^2;
basinPop_ppl=sum(datain.Pop(bidx))*(cellsz_m/1e3)^2 ; 
if plotallfigs
    %% subcatchment map
    figure
    subplot(2,2,1:2)
    h=imagescnan(maskBasin(catchments,catchments>0&catchments<109));
    colormap(cmap8)
    colorbar%'Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'Ticklabels',(basindata.basinnames),'TickDirection','none')    % savefig(fullfile(sprintf("%s/%s/TheoreticalPotential_SubbasinDefintion.fig",root,res)))
    %set(h, 'AlphaData',0.9)
    subplot(2,2,3)
    bar(diag(basindata.Area_m2),'stacked')
    ylabel("Area (m^2)")
    applymyplotformat('',cmap8)
    xticklabels(basindata.basinnames)
    
    subplot(2,2,4)
    bar(diag(basindata.Pop_ppl),'stacked')
    ylabel("# of inhabitants")
    applymyplotformat('',cmap8)
    xticklabels(basindata.basinnames)
end
disp("Loaded basin info for Basin")

%% Get subbasin-wise VISUALIZED potential
load(fullfile(rootf,'data\UI\data\Existing+UC+All_3.mat'), 'existing_dams','existing_reservoirs')
viz_status=existing_dams.status;
viz_GWh=existing_dams.GWh;
viz_c=existing_dams.c_dams;
viz_r=existing_dams.r_dams;

isel=sub2ind(size(datain.Regions), viz_r, viz_c);
viz_subIDs=datain.Regions(isel);

for st=basindata.basinIDs'
    selidxs=viz_subIDs ==st;
    viz_subPot_GWh(st-100,:) = sum(viz_GWh(selidxs));
    viz_subPot_num(st-100,:)=sum(selidxs);
    
    % sum of projects in existing state
    viz_subPot_GWh_exist(st-100,:) = sum(viz_GWh(selidxs & strcmp(viz_status,'Existing')));
    viz_subPot_num_exist(st-100,:) = sum(selidxs & strcmp(viz_status,'Existing'));
    
    % sum of projects in not existing state
    viz_subPot_GWh_remain(st-100,:) = sum(viz_GWh(selidxs & ~strcmp(viz_status,'Existing')));
end

viz_compiled=existing_dams(:,1:7);
viz_compiled.subbasinIDs=viz_subIDs;
viz_compiled.subbasins= categorical(viz_subIDs,basindata.basinIDs,basindata.basinnames);

% Add basin wide total to end of subPot
viz_subPot_GWh(9,:) = sum(viz_GWh);
viz_subPot_num(9,:)=numel(viz_GWh);

viz_subPot_GWh_exist(9,:) = sum(viz_GWh( strcmp(viz_status,'Existing')));
viz_subPot_num_exist(9,:) = sum(strcmp(viz_status,'Existing'));
viz_subPot_GWh_remain(9,:) = sum(viz_GWh(~strcmp(viz_status,'Existing')));

%% Get subbasin-wise THEORETICAL potential
load(fullfile(rootof,'\fig_theoretical_fixedDEM_QbasedChannel\500m\TheoreticalPot.mat'), 'channel_basinEnergy')
ch_basinEnergy=channel_basinEnergy(:,:,1);
for st=basindata.basinIDs'
    selidxs=ch_basinEnergy > 0 & catchments ==st;
    theory_subPot_GWh(st-100,:) = sum(ch_basinEnergy(selidxs),'all');
    theory_subPot_num(st-100,:)=sum(selidxs,'all');
end

% Add basin wide total to end of subPot
selidxs=ch_basinEnergy > 0 ;
theory_subPot_GWh(9,:) = sum(ch_basinEnergy(selidxs),'all');
theory_subPot_num(9,:)=sum(selidxs,'all');

%Plot theoretical and visualized potential
figure;
bar(theory_subPot_GWh,'FaceAlpha',0.5)
hold all
bar(viz_subPot_GWh,'FaceAlpha',0,'EdgeColor','r')
bar(viz_subPot_GWh_exist,'FaceAlpha',0,'EdgeColor','g')
legend(pottypes6([1,5,6]))
    xticklabels([basindata.basinnames{:},"All Indus"])
%% Loop through and compile outputs for FULL 3 policy types
runcnt=1;
runname_prefix='R103_Energy';
for scen={'Full', 'Remain'}%   'Full';  % 'Remain' ;  %
    for policyname={'Large','Medium','Mixed'}
        for consttype={'Tech_Fin', 'Sust_RiskAverse'}
            % Separate tech and econ data
            runnames{runcnt}=strjoin([{runname_prefix},scen(:)', policyname(:)', consttype(:)'],'_');
            % Prepare exact filename for matfile
            pcsout{runcnt} = load(fullfile(rootf, 'output',scen{:}, continent_in,runnames{runcnt},...
                sprintf('COEPOT_b%d_do.mat', nbasin)));
            pcsout{runcnt}.runname=runnames{runcnt};
            
            pcsout{runcnt}=rmfield(pcsout{runcnt}, {'RDDepth','RDCountry_id','RDRegion_id'} );
            runcnt=runcnt+1;
        end
    end
end
disp("Loaded output data for Basin - Energy ")

%% Loop through and compile outputs for FULL Sust Hazard Rep scenarios
runcnthaz=1;
runname_prefix='R103_HazRep';
for scen={'Full'}%   'Full';  % 'Remain' ;  %
    for policyname={'Mixed'}
        for consttype={'Sust_Costbased', 'Sust_Multi-hazard'}
            % Separate tech and econ data
            runnameshaz{runcnthaz}=strjoin([{runname_prefix},scen(:)', policyname(:)', consttype(:)'],'_');
            % Prepare exact filename for matfile
            pcsouthaz{runcnthaz} = load(fullfile(rootf, 'output',scen{:}, continent_in,runnameshaz{runcnthaz},...
                sprintf('COEPOT_b%d_do.mat', nbasin)));
            pcsouthaz{runcnthaz}.runname=runnameshaz{runcnthaz};
            pcsouthaz{runcnthaz}=rmfield(pcsouthaz{runcnthaz}, {'RDDepth','RDCountry_id','RDRegion_id'} );
            runcnthaz=runcnthaz+1;
        end
    end
end
disp("Loaded output data for Basin - Hazard Rep")

%% Loop through and compile outputs for FULL Tech2sust scenarios
runcntT2S=1;
runname_prefix='R103_Tech2Sust';
consttype={'Tech_Fin'};
for scen={'Full'}%   'Full';  % 'Remain' ;  %
    for policyname={'Mixed'}
        for addconst={'+wc','+wc+eflow','+wc+eflow+PA'}
            % Separate tech and econ data
            runnamesT2S{runcntT2S}=strjoin([[runname_prefix, addconst{:}],scen(:)', policyname(:) consttype(:)'],'_');
            % Prepare exact filename for matfile
            pcsoutT2S{runcntT2S} = load(fullfile(rootf, 'output',scen{:}, continent_in,runnamesT2S{runcntT2S},...
                sprintf('COEPOT_b%d_do.mat', nbasin)));
            pcsoutT2S{runcntT2S}.runname=runnamesT2S{runcntT2S};
            pcsoutT2S{runcntT2S}=rmfield(pcsoutT2S{runcntT2S}, {'RDDepth','RDCountry_id','RDRegion_id'} );
            runcntT2S=runcntT2S+1;
        end
    end
end
disp("Loaded output data for Basin - Tech2Sust")

%% Assign Sub-basin and country to each project
% ENERGY
for i=1:numel(pcsout)
    % odd i vals are tech/econ, rest are sust
    mydata=pcsout{i};
    mydata.subbasinIDs=nan(size(mydata.ross));
    mydata.countryIDs=nan(size(mydata.ross));
    
    % Correct data only for not nan cases
    valid_idxs=find(~isnan(mydata.ross));
    rsel = mydata.ross(valid_idxs);
    csel = mydata.coss(valid_idxs);
    
    isel=sub2ind(size(datain.Regions), rsel, csel);
    mydata.subbasinIDs(valid_idxs)=datain.Regions(isel);
    mydata.countryIDs(valid_idxs)=datain.Countries(isel);
    pcsout{i}=mydata;
end

% HAZARD
for i=1:numel(pcsouthaz)
    % odd i vals are tech/econ, rest are sust
    mydata=pcsouthaz{i};
    mydata.subbasinIDs=nan(size(mydata.ross));
    mydata.countryIDs=nan(size(mydata.ross));
    
    % Correct data only for not nan cases
    valid_idxs=find(~isnan(mydata.ross));
    rsel = mydata.ross(valid_idxs);
    csel = mydata.coss(valid_idxs);
    
    isel=sub2ind(size(datain.Regions), rsel, csel);
    mydata.subbasinIDs(valid_idxs)=datain.Regions(isel);
    mydata.countryIDs(valid_idxs)=datain.Countries(isel);
    pcsouthaz{i}=mydata;
end

%% Assign Sub-basin and take sub-basin sums for Tech2Sust
for runid=1:numel(pcsoutT2S)
    pcsoutT2S{runid}.subbasinIDs=getRegionID(pcsoutT2S{runid}.ross, pcsoutT2S{runid}.coss, datain.Regions);
    pcsoutT2S{runid}.countryIDs=getRegionID(pcsoutT2S{runid}.ross, pcsoutT2S{runid}.coss, datain.Countries);

    valid_idxs=~isnan(pcsoutT2S{runid}.ross); % This is same for coss
    % Take subbasin sums
    % this is for tech pot but tech2 sust is only after financial filters so skip
    %subPot_GWh_T2S(:,runid)=getRegionalSum(101:108,pcsoutT2S{runid}.subbasinIDs(valid_idxs), pcsoutT2S{runid}.PnetAlls(valid_idxs)); 
    selecon= valid_idxs & pcsoutT2S{runid}.COEAlls<=costlim;
    subPot_GWh_tmp(:,runid)=getRegionalSum(101:108,pcsoutT2S{runid}.subbasinIDs(selecon), pcsoutT2S{runid}.PnetAlls(selecon));
end
subPot_GWh_T2S=array2table(subPot_GWh_tmp,'VariableNames',strcat('Tech+Fin',{'+wc','+wc+eflow','+wc+eflow+PA'}'));
subPot_GWh_T2S.basinID=basinlabels.basinIDs(1:8);
subPot_GWh_T2S.basinname=basinlabels.basinnames(1:8);
disp("Assigned subbasin and country to each point - Tech2Sust")

%Save for CAS
ofname=fullfile(rootof,'Figs_trial','SubBasinData_Tech2Sust.xlsx');
writetable(subPot_GWh_T2S(:,["basinID","basinname","Tech+Fin+wc+eflow","Tech+Fin+wc+eflow+PA"]),ofname,'Sheet','GWhperYr')

%% Calculate Sub-basin wise potential sums including existing - for main scenarios
sc_idx=0;  % scenarios are layered in the 3Dmatrix
for runid=1:length(runnames)
    % Load and process techecon data
    mydata=pcsout{runid};
    if rem(runid,2)>0 % start new matrix sheet for odd runs and add theorypot to it
        sc_idx=sc_idx+1;
        % Add theory pot to start
        subPot_GWh(:,1,sc_idx)=theory_subPot_GWh;
        subPot_num(:,1,sc_idx)=theory_subPot_num;
        % Add vis pot to end
        subPot_GWh(:,5,sc_idx)=viz_subPot_GWh; % Full
        subPot_num(:,5,sc_idx)=viz_subPot_num;
        % Add exist pot to end
        subPot_GWh(:,6,sc_idx)=viz_subPot_GWh_exist; % Full
        subPot_num(:,6,sc_idx)=viz_subPot_num_exist;
        
        if sc_idx>3 %For remaining case only add vis pot, theory pot is not relevant
            subPot_GWh(:,5,sc_idx)=viz_subPot_GWh_remain; % Remain
        end
    end
    
    % get not nan datas
    valid_idxs=find(~isnan(mydata.ross)); % This is same for cos_rdid
    GWhsel= mydata.PnetAlls(valid_idxs);
    COEsel= mydata.COEAlls(valid_idxs);
    subIDs=mydata.subbasinIDs(valid_idxs);
    
    % Get subbasin-wise tech econ sust potential
    for st=101:108
        if rem(runid,2)>0 % for odd sc save both tech and econ in current row
            % Tech pot
            selidxs=subIDs ==st;
            subPot_GWh(st-100,2,sc_idx) = sum(GWhsel(selidxs));
            subPot_num(st-100,2,sc_idx)=sum(selidxs);
            % Econ pot
            selidxs=subIDs ==st & COEsel<=costlim;
            subPot_GWh(st-100,3,sc_idx) = sum(GWhsel(selidxs) );
            subPot_num(st-100,3,sc_idx) = sum(selidxs);
        else  % for even sc save sust in previous row
            % Sust pot
            selidxs=subIDs ==st & COEsel<=costlim;
            subPot_GWh(st-100,4,sc_idx) = sum(GWhsel(selidxs));
            subPot_num(st-100,4,sc_idx)=sum(selidxs);
        end
    end
    
    % Get basin-wide tech econ sust potential
    if rem(runid,2)>0 % for odd sc save both tech and econ in current row
        % Tech pot
        selidxs=subIDs>0;
        subPot_GWh(end,2,sc_idx) = sum(GWhsel(selidxs));
        subPot_num(end,2,sc_idx)=sum(selidxs);
        % Econ pot
        selidxs= COEsel<=costlim;
        subPot_GWh(end,3,sc_idx) = sum(GWhsel(selidxs) );
        subPot_num(end,3,sc_idx) = sum(selidxs);
    else  % for even sc save sust in previous row
        % Sust pot
        selidxs= COEsel<=costlim;
        subPot_GWh(end,4,sc_idx) = sum(GWhsel(selidxs));
        subPot_num(end,4,sc_idx)=sum(selidxs);
    end
end

disp("Calculated sub-basin wise totals for main scenarios")

%% Get energy per area and per capita for subbasin
for sc_idx=1:size(subPot_num,3)
    % For all subbasins
    for st=101:108
        subPot_GWh_perarea(st-100,:,sc_idx) =  subPot_GWh(st-100,:,sc_idx) /basindata.Area_m2(basindata.basinIDs==st);
        subPot_GWh_percapita(st-100,:,sc_idx) = subPot_GWh(st-100,:,sc_idx) /basindata.Pop_ppl(basindata.basinIDs==st);
    end
    % For all basin
subPot_GWh_perarea(9,:,sc_idx) =  subPot_GWh(9,:,sc_idx)/basinArea_m2;
subPot_GWh_percapita(9,:,sc_idx) = subPot_GWh(9,:,sc_idx)/basinPop_ppl;
end
disp("Calculated per area and per capita vals")


%% Plot Sub-basin wise potential sums
if plotallfigs
%% Plot all 6 pot types simply
figure
selpots=1:size(subPot_GWh,2);
selbasins=1:9;
cmap=lines(6);%magma(8);%cbrewer2('PuOR',length(selcols)+1); % Colorblind + printer friendly %lines;
cmap3=cmap(selpots,:);
for sc_idx=1:3
    subplot(3,3,sc_idx)
    bar(subPot_GWh(selbasins,selpots,sc_idx)/1000,'FaceAlpha',baralpha,'EdgeAlpha',baralpha)
    if sc_idx==1;   ylabel({'Total potential';'(TWh/yr)'},'FontWeight','bold');    end
    xticklabels(label_subbasin_basin)
    xtickangle(30)
    applymyplotformat(searchtypes{sc_idx},cmap3)
    set(gca,'XGrid','off')
    %ylim([0, 72])
    
    subplot(3,3,sc_idx+3)
    bar(subPot_GWh_perarea(selbasins,selpots,sc_idx)*1e6,'FaceAlpha',baralpha,'EdgeAlpha',baralpha)
    if sc_idx==1;   ylabel({'Specific potential'; '(kWh/yr per m^2)'},'FontWeight','bold'); end
    xticklabels(label_subbasin_basin)
    xtickangle(30)
    applymyplotformat('',cmap3)
    set(gca,'XGrid','off')
    %ylim([0, 2.6])
    
    subplot(3,3,sc_idx+6)
    bar(subPot_GWh_percapita(selbasins,selpots,sc_idx)*1e3,'FaceAlpha',baralpha,'EdgeAlpha',baralpha)
    if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
    xtickangle(30)
    xticklabels(label_subbasin_basin)
    applymyplotformat('',cmap3)
    set(gca,'XGrid','off')
    %ylim([0, 35])    
end
subplot(3,3,1)
legend(pottypes6(selpots),'Orientation','horizontal','Location','north') %,'NumColumns',4  )

    %% Plot Sub-basin totals - all Tech/econ/sust/vis/exist pot types as bars, add theory pot as empty bars
figure
selbasins=1:9;
selpots=2:size(subPot_GWh,2);
cmap=lines(7);%magma(8);%cbrewer2('PuOR',length(selcols)+1); % Colorblind + printer friendly %lines;
cmap3=cmap(selpots,:);
for sc_idx=1:3
    subplot(1,3,sc_idx)
    bar(subPot_GWh(selbasins,selpots,sc_idx)/1000,'FaceAlpha',baralpha,'EdgeAlpha',0)
    if sc_idx==1;   ylabel({'Total potential';'(TWh/yr)'},'FontWeight','bold');    end
    xticklabels(label_subbasin_basin)
    xtickangle(30)
    applymyplotformat(searchtypes{sc_idx},cmap3)
    set(gca,'XGrid','off')
    %ylim([0, 72])
    hold all
    % add theoretical pot on top
    b12=bar(subPot_GWh(selbasins,1,sc_idx)/1000,'FaceColor',cmap(1,:),'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'LineStyle','-');
    set(gca, 'YScale', 'log')    
end
legend(pottypes6{[selpots 1]},'Orientation','horizontal','Location','north') %,'NumColumns',4  )

%% Plot all Tech/econ/sust pot types as bars, add basinwide and existing as empty bars
figure
selpots=2:4; %size(subPot_GWh,2);
cmap3=cmap(selpots,:);
for sc_idx=1:3
    subplot(2,3,sc_idx)
    idata=subPot_GWh_perarea*1e6; % subPot_GWh_percapita*1e3; % 
    bar(idata(selbasins,selpots,sc_idx),'FaceAlpha',baralpha,'EdgeAlpha',0)
    if sc_idx==1;   ylabel({'Specific potential'; '(kWh/yr per m^2)'},'FontWeight','bold'); end
    xtickangle(30)
    xticklabels(label_subbasin_basin)
    applymyplotformat('',cmap3)
    set(gca,'XGrid','off')
    %ylim([0, 35]) 
    
    hold on
    b2=bar(idata(selbasins,5,sc_idx),'LineStyle','-','FaceColor',[1 1 1]*.8,'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'DisplayName','Visualized');
    b3=bar(idata(selbasins,6,sc_idx),'LineStyle','-','FaceColor','r','EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'DisplayName','Existing');
    %b4=bar(repmat(idata(9,selpots,sc_idx),8,1),'LineStyle','-','FaceColor','k','EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'DisplayName','Basin level value');
    
    subplot(2,3,sc_idx+3)
    idata=subPot_GWh_percapita*1e3; % subPot_GWh_percapita*1e3; % 
    bar(idata(selbasins,selpots,sc_idx),'FaceAlpha',baralpha,'EdgeAlpha',0)
    if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
    xtickangle(30)
    xticklabels(label_subbasin_basin)
    applymyplotformat('',cmap3)
    set(gca,'XGrid','off')
    %ylim([0, 35]) 
    
    hold on
    b21=bar(idata(selbasins,5,sc_idx),'LineStyle','-','FaceColor',[1 1 1]*.8,'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'DisplayName','Visualized');
    b31=bar(idata(selbasins,6,sc_idx),'LineStyle','-','FaceColor','r','EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'DisplayName','Existing');
    %b41=bar(repmat(idata(9,selpots,sc_idx),8,1),'LineStyle','-','FaceColor','k','EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'DisplayName','Basin level value');
        
   xtickangle(30)
    xticklabels(label_subbasin_basin)
    applymyplotformat('',cmap3)
    set(gca,'XGrid','off')
    %ylim([0, 35])    
end
%    legend(pottypes6(selpots),'Orientation','horizontal','Location','north') %,'NumColumns',4  )
legend([pottypes6{[selpots 5 6]}, "All basin value"],'Orientation','horizontal','Location','north') %,'NumColumns',4  )

%% Onescenario: Bar plot w existing and viz pot one on top of another
figure
selpots=1:4;%size(subPot_GWh,2);
 bar(subPot_GWh(selbasins,selpots,sc_idx)/1000,'FaceAlpha',baralpha,'EdgeAlpha',0)
    if sc_idx==1;   ylabel({'Total potential';'(TWh/yr)'},'FontWeight','bold');    end
    xticklabels(label_subbasin_basin)
    xtickangle(30)
    applymyplotformat(searchtypes{sc_idx},cmap3)
    set(gca,'XGrid','off')
    %ylim([0, 72])
    hold all
    % add existing and viz pot on top
    b2=bar(subPot_GWh(selbasins,5,sc_idx)/1000,'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'LineStyle','-');
    b3=bar(subPot_GWh(selbasins,6,sc_idx)/1000,'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'LineStyle','-');
    set(gca, 'YScale', 'log')
    legend(pottypes6([selpots 5 6]),'Orientation','horizontal','Location','north') %,'NumColumns',4  )

%% Onescenario: Plot all basin totals and existing,visualized val over the others
selbasins=1:8;
selpots=1:4;
idata=subPot_GWh_percapita*1e3; % 
figure
    bar(idata(selbasins,selpots,sc_idx),'FaceAlpha',baralpha,'EdgeAlpha',0)
    if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
    xtickangle(30)
    xticklabels(label_subbasin_basin)
    applymyplotformat('',cmap3)
    set(gca,'XGrid','off')
    %ylim([0, 35]) 
    
    hold on
    b2=bar(idata(selbasins,5,sc_idx),'LineStyle','-','EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'DisplayName','Visualized');
    b3=bar(idata(selbasins,6,sc_idx),'LineStyle','-','FaceColor','r','EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'DisplayName','Existing');
    b4=bar(repmat(idata(9,selpots,sc_idx),8,1),'LineStyle','--','FaceColor','k','EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'DisplayName','Basin level value');
    legend([pottypes6{[selpots 5 6]}, "All basin value"],'Orientation','horizontal','Location','north') %,'NumColumns',4  )


    %% Stacked basins bar of FULL potential scenarios - only TFS
    scnames=strrep(extractBetween(runnames([1:2:12]),'Energy_','_Tech'),'_', '-');
    figure
    selpots=2:4;
    for sc_idx=1:3
        subplot(3,3,sc_idx)
        bar(subPot_GWh(:,selpots,sc_idx)'/1000,'stacked')
        if sc_idx==1;   ylabel({'Total potential';'(TWh/yr)'},'FontWeight','bold');    end
        %xtickangle(45)
        xticklabels([])
        applymyplotformat(scnames{sc_idx},cmap8)
        ylim([0, 310])
        
        subplot(3,3,sc_idx+3)
        bar(subPot_GWh_perarea(:,selpots,sc_idx)'*1e6,'stacked')
        if sc_idx==1;   ylabel({'Specific potential'; '(kWh/yr per m^2)'},'FontWeight','bold'); end
        xticklabels([])
        %xtickangle(45)
        applymyplotformat('',cmap8)
        ylim([0, 12])
        
        subplot(3,3,sc_idx+6)
        h=bar(subPot_GWh_percapita(:,selpots,sc_idx)'*1e3,'stacked');
        if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
        xticklabels(pottypes5_short(selpots))
        %xtickangle(45)
        applymyplotformat('',cmap8)
        ylim([0, 100])
        
    end
    %legend(flip(h), flip(basindata.basinnames),'Location','north' ) %,'NumColumns',4
    legend(flip(h), flip(basindata.basinnames),'Orientation','horizontal','Location','north','Position',[0.2 0.9 0.5 0.03], 'NumColumns',4)
    
    %% Not stacked potential type bar of FULL potential scenarios  - only TFS
    figure
    selpots=2:4; % - only TFS
    cmap3=cbrewer2('PuOR',length(selpots)+1); % Colorblind + printer friendly %lines;
    cmap3=cmap3([1,3,4],:);
    for sc_idx=1:3
        subplot(3,3,sc_idx)
        bar(subPot_GWh(:,selpots,sc_idx)/1000)
        if sc_idx==1;   ylabel({'Total potential';'(TWh/yr)'},'FontWeight','bold');    end
        xticklabels([])
        %xtickangle(45)
        applymyplotformat(scnames{sc_idx},cmap3)
        ylim([0, 72])
        
        subplot(3,3,sc_idx+3)
        bar(subPot_GWh_perarea(:,selpots,sc_idx)*1e6)
        if sc_idx==1;   ylabel({'Specific potential'; '(kWh/yr per m^2)'},'FontWeight','bold'); end
        xticklabels([])
        %xtickangle(45)
        applymyplotformat('',cmap3)
        ylim([0, 2.6])
        
        subplot(3,3,sc_idx+6)
        bar(subPot_GWh_percapita(:,selpots,sc_idx)*1e3)
        if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
        %xtickangle(45)
        applymyplotformat('',cmap3)
        xticklabels(basindata.basinnames)
        ylim([0, 35])
        
    end
    % legend(basindata.basinnames,'Orientation','horizontal','Location','north','NumColumns',4  )
    legend(pottypes5(selpots),'Orientation','horizontal','Location','north') %,'NumColumns',4  )
    %% Stacked basins bar of FULL potential scenarios - theory + TFS
    figure
    selpots=1:4;
    scnames=strrep(extractBetween(runnames([1:2:12]),'Energy_','_Tech'),'_', '-');
    for sc_idx=1:3
        subplot(3,3,sc_idx)
        bar(subPot_GWh(:,selpots,sc_idx)'/1000,'stacked')
        if sc_idx==1;   ylabel({'Total potential';'(TWh/yr)'},'FontWeight','bold');    end
        xticklabels([])
        applymyplotformat(scnames{sc_idx},cmap8)
        %ylim([0, 300])
        %xtickangle(45)
        
        subplot(3,3,sc_idx+3)
        bar(subPot_GWh_perarea(:,selpots,sc_idx)'*1e6,'stacked')
        if sc_idx==1;   ylabel({'Specific potential'; '(kWh/yr per m^2)'},'FontWeight','bold'); end
        xticklabels([])
        applymyplotformat('',cmap8)
        %ylim([0, 12])
        
        subplot(3,3,sc_idx+6)
        h=bar(subPot_GWh_percapita(:,selpots,sc_idx)'*1e3,'stacked');
        if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
        xticklabels(pottypes5_short(selpots))
        %xtickangle(45)
        applymyplotformat('',cmap8)
        %ylim([0, 100])
        
    end
    %legend(flip(h), flip(basindata.basinnames),'Location','north' ) %,'NumColumns',4
    legend(flip(h), flip(basindata.basinnames),'Orientation','horizontal','Location','north','Position',[0.2 0.9 0.5 0.03], 'NumColumns',4)
    
    
    %% Not stacked potential type bar of FULL potential scenarios- theory + TFS
    selpots=1:4; % theory + TFS
    figure
    cmap3=cbrewer2('PuOR',length(selpots)); % Colorblind + printer friendly %lines;
    %cmap3=cmap3(2:end,:);
    for sc_idx=1:3
        subplot(3,3,sc_idx)
        bar(subPot_GWh(:,selpots,sc_idx)/1000)
        if sc_idx==1;   ylabel({'Total potential';'(TWh/yr)'},'FontWeight','bold');    end
        xticklabels([])
        %xtickangle(45)
        applymyplotformat(scnames{sc_idx},cmap3)
        %ylim([0, 70])
        
        subplot(3,3,sc_idx+3)
        bar(subPot_GWh_perarea(:,selpots,sc_idx)*1e6)
        if sc_idx==1;   ylabel({'Specific potential'; '(kWh/yr per m^2)'},'FontWeight','bold'); end
        xticklabels([])
        %xtickangle(45)
        applymyplotformat('',cmap3)
        %ylim([0, 2.52])
        
        subplot(3,3,sc_idx+6)
        bar(subPot_GWh_percapita(:,selpots,sc_idx)*1e3)
        if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
        %  xticklabels(["Technical", "Financial", "Sustainable"])
        %xtickangle(45)
        applymyplotformat('',cmap3)
        xticklabels(basindata.basinnames)
        %ylim([0, 35])
    end
    % legend(basindata.basinnames,'Orientation','horizontal','Location','north','NumColumns',4  )
    legend(pottypes5(selpots),'Orientation','horizontal','Location','north') %,'NumColumns',4  )
    
    %% Stacked basins bar of FULL potential scenarios - theory + TFS + Visualized
    figure
    selpots=1:5; %- theory + TFS
    for sc_idx=1:3
        subplot(3,3,sc_idx)
        bar(subPot_GWh(:,selpots,sc_idx)'/1000,'stacked')
        if sc_idx==1;   ylabel({'Total potential';'(TWh/yr)'},'FontWeight','bold');    end
        xticklabels([])
        applymyplotformat(scnames{sc_idx},cmap8)
        %ylim([0, 300])
        %xtickangle(45)
        
        subplot(3,3,sc_idx+3)
        bar(subPot_GWh_perarea(:,selpots,sc_idx)'*1e6,'stacked')
        if sc_idx==1;   ylabel({'Specific potential'; '(kWh/yr per m^2)'},'FontWeight','bold'); end
        xticklabels([])
        applymyplotformat('',cmap8)
        %ylim([0, 12])
        
        subplot(3,3,sc_idx+6)
        h=bar(subPot_GWh_percapita(:,selpots,sc_idx)'*1e3,'stacked');
        if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
        xticklabels(pottypes5_short(selpots))
        %xtickangle(45)
        applymyplotformat('',cmap8)
        %ylim([0, 100])
        
    end
    %legend(flip(h), flip(basindata.basinnames),'Location','north' ) %,'NumColumns',4
    legend(flip(h), flip(basindata.basinnames),'Orientation','horizontal','Location','north','Position',[0.2 0.9 0.5 0.03], 'NumColumns',4)
    
    
    %% Not stacked potential type bar of FULL potential scenarios - theory + TFS
    selpots=1:5; %- theory + TFS
    figure
    cmap3=cbrewer2('PuOR',length(selpots)); % Colorblind + printer friendly %lines;
    %cmap3=cmap3(2:end,:);
    for sc_idx=1:3
        subplot(3,3,sc_idx)
        bar(subPot_GWh(:,selpots,sc_idx)/1000)
        if sc_idx==1;   ylabel({'Total potential';'(TWh/yr)'},'FontWeight','bold');    end
        xticklabels([])
        %xtickangle(45)
        applymyplotformat(scnames{sc_idx},cmap3)
        %ylim([0, 70])
        
        subplot(3,3,sc_idx+3)
        bar(subPot_GWh_perarea(:,selpots,sc_idx)*1e6)
        if sc_idx==1;   ylabel({'Specific potential'; '(kWh/yr per m^2)'},'FontWeight','bold'); end
        xticklabels([])
        %xtickangle(45)
        applymyplotformat('',cmap3)
        %ylim([0, 2.52])
        
        subplot(3,3,sc_idx+6)
        bar(subPot_GWh_percapita(:,selpots,sc_idx)*1e3)
        if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
        %  xticklabels(["Technical", "Financial", "Sustainable"])
        %xtickangle(45)
        applymyplotformat('',cmap3)
        xticklabels(basindata.basinnames)
        %ylim([0, 35])
    end
    % legend(basindata.basinnames,'Orientation','horizontal','Location','north','NumColumns',4  )
    legend(pottypes5(selpots),'Orientation','horizontal','Location','north') %,'NumColumns',4  )
    
    %% Stacked basins bar of FULL potential scenarios -  TFS + Visualized
    figure
    selpots=2:5; %-  TFS + Visualized
    for sc_idx=1:3
        subplot(3,3,sc_idx)
        bar(subPot_GWh(:,selpots,sc_idx)'/1000,'stacked')
        if sc_idx==1;   ylabel({'Total potential';'(TWh/yr)'},'FontWeight','bold');    end
        xticklabels([])
        applymyplotformat(scnames{sc_idx},cmap8)
        ylim([0, 310])
        %xtickangle(45)
        
        subplot(3,3,sc_idx+3)
        bar(subPot_GWh_perarea(:,selpots,sc_idx)'*1e6,'stacked')
        if sc_idx==1;   ylabel({'Specific potential'; '(kWh/yr per m^2)'},'FontWeight','bold'); end
        xticklabels([])
        applymyplotformat('',cmap8)
        ylim([0, 11])
        
        subplot(3,3,sc_idx+6)
        h=bar(subPot_GWh_percapita(:,selpots,sc_idx)'*1e3,'stacked');
        if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
        xticklabels(pottypes5_short(selpots))
        %xtickangle(45)
        applymyplotformat('',cmap8)
        ylim([0, 90])
        
    end
    %legend(flip(h), flip(basindata.basinnames),'Location','north' ) %,'NumColumns',4
    legend(flip(h), flip(basindata.basinnames),'Orientation','horizontal','Location','north','Position',[0.2 0.9 0.5 0.03], 'NumColumns',4)
        
    %% Not stacked potential type bar of FULL potential scenarios -  TFS + Visualized
    selpots=2:5; %-  TFS + Visualized
    
    figure
    cmap3=cbrewer2('PuOR',length(selpots)); % Colorblind + printer friendly %lines;
    %cmap3=cmap3(2:end,:);
    for sc_idx=1:3
        subplot(3,3,sc_idx)
        bar(subPot_GWh(:,selpots,sc_idx)/1000)
        if sc_idx==1;   ylabel({'Total potential';'(TWh/yr)'},'FontWeight','bold');    end
        xticklabels([])
        %xtickangle(45)
        applymyplotformat(scnames{sc_idx},cmap3)
        ylim([0, 165])
        
        subplot(3,3,sc_idx+3)
        bar(subPot_GWh_perarea(:,selpots,sc_idx)*1e6)
        if sc_idx==1;   ylabel({'Specific potential'; '(kWh/yr per m^2)'},'FontWeight','bold'); end
        xticklabels([])
        %xtickangle(45)
        applymyplotformat('',cmap3)
        ylim([0, 2.7])
        
        subplot(3,3,sc_idx+6)
        bar(subPot_GWh_percapita(:,selpots,sc_idx)*1e3)
        if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
        %  xticklabels(["Technical", "Financial", "Sustainable"])
        %xtickangle(45)
        applymyplotformat('',cmap3)
        xticklabels(basindata.basinnames)
        ylim([0, 33])
    end
    % legend(basindata.basinnames,'Orientation','horizontal','Location','north','NumColumns',4  )
    legend(pottypes5(selpots),'Orientation','horizontal','Location','north') %,'NumColumns',4  )
end

%% Calculate Sub-basin wise potential sums separated by RP and DP only for Mixed Full+Remaining case for plotting in R
sc_idx=0;
%SysIDs 1=DP, 2=RP
for runid=1:length(runnames)
    % Load and process techecon data
    mydata=pcsout{runid};
    if rem(runid,2)>0 % start new matrix sheet for odd runs and add theory/vispot to it
        sc_idx=sc_idx+1; %increment row for Tech/Econ scenario, not for Sust scenario
        % Add theory pot to start
        subPot_GWh_RP(:,1,sc_idx)=theory_subPot_GWh(1:8);
        subPot_GWh_DP(:,1,sc_idx)=nan;
        
        % Add vis pot to end of RP data
        subPot_GWh_RP(:,5,sc_idx)=viz_subPot_GWh(1:8);
        subPot_GWh_DP(:,5,sc_idx)=nan;
        
        % Add vis pot to end of RP data
        subPot_GWh_RP(:,5,sc_idx)=viz_subPot_GWh(1:8);
        subPot_GWh_DP(:,5,sc_idx)=nan;
        
        if sc_idx>3
            subPot_GWh_RP(:,5,sc_idx)=viz_subPot_GWh_remain(1:8);
            subPot_GWh_DP(:,5,sc_idx)=nan;
            %plant type not valid for vis data
        end
    end
    
    % get not nan datas
    valid_idxs=find(~isnan(mydata.ross)); % This is same for cos_rdid
    GWhsel= mydata.PnetAlls(valid_idxs);
    COEsel= mydata.COEAlls(valid_idxs);
    subIDs=mydata.subbasinIDs(valid_idxs);
    plantID=mydata.SysIDs(valid_idxs);
    
    % Get subbasin-wise tech econ sust potential
    for st=101:108
        if rem(runid,2)>0 % for odd runid save both tech and econ
            % Tech pot
            selidxs=subIDs ==st;
            subPot_GWh_RP(st-100,2,sc_idx) = sum(GWhsel(selidxs & plantID ==2));
            subPot_GWh_DP(st-100,2,sc_idx) = sum(GWhsel(selidxs & plantID ==1));
                        % Econ pot
            selidxs=subIDs ==st & COEsel<=costlim;
            subPot_GWh_RP(st-100,3,sc_idx) = sum(GWhsel(selidxs& plantID ==2) );
            subPot_GWh_DP(st-100,3,sc_idx) = sum(GWhsel(selidxs& plantID ==1) );
        else  % for even runid save only sust
            % sust pot
            selidxs=subIDs ==st & COEsel<=costlim;
            subPot_GWh_RP(st-100,4,sc_idx) = sum(GWhsel(selidxs& plantID ==2));
            subPot_GWh_DP(st-100,4,sc_idx) = sum(GWhsel(selidxs& plantID ==1));
        end
    end
end

%% Prep data for plotting in R
subPot_long=table;
scnames=strrep(extractBetween(runnames([1:2:12]),'Energy_','_Tech'),'_', '-');

for sc_idx=3 %1:3  % Run only remain
    %get full scenario data
    RPdata=subPot_GWh_RP(:,:,sc_idx);  % basin (8) x potential type (5) x scen (Full then Remain)
    DPdata=subPot_GWh_DP(:,:,sc_idx);
    % create new table and concatenate
    dcompile=array2table([RPdata(:); DPdata(:)],'VariableNames', "Full_GWh_year");
    dcompile.planttype = repelem(planttypes,numel(RPdata(:)))';
    dcompile.pottype = categorical(repmat(repelem(pottypes5,length(basindata.basinnames))',2,1),pottypes5([5 1:4]));
    dcompile.basin = categorical(repmat(basindata.basinnames, size(RPdata,2)*2,1),basindata.basinnames);
    dcompile.scname = repmat(scnames(sc_idx),height(dcompile),1);
    % get remain scenario data
    RPdata=subPot_GWh_RP(:,:,sc_idx+3);  % basin (8) x potential type (5) x scen (Full then Remain)
    DPdata=subPot_GWh_DP(:,:,sc_idx+3);
    
    dcompile=addvars(dcompile(:,[2:4 5 1]),[RPdata(:); DPdata(:)],'NewVariableNames', "Remain_GWh_year");
    
    subPot_long=[subPot_long; dcompile];
end

%% Save for Saurav
if save4saurav
    for ssel=1:3
    scfname=fullfile(rootof,'Figs_trial',sprintf('SubBasinData_%s.xlsx',searchtypes{ssel}));
    tbl_GWh=table(subPot_GWh(:,1,ssel),subPot_GWh(:,2,ssel), subPot_GWh(:,3,ssel),...
        subPot_GWh(:,4,ssel), subPot_GWh(:,5,ssel), subPot_GWh(:,6,ssel),'VariableNames',pottypes6,'RowNames',label_subbasin_basin);
    
    tbl_perarea=table(subPot_GWh_perarea(:,1,ssel),subPot_GWh_perarea(:,2,ssel), subPot_GWh_perarea(:,3,ssel),...
        subPot_GWh_perarea(:,4,ssel), subPot_GWh_perarea(:,5,ssel), subPot_GWh_perarea(:,6,ssel),'VariableNames',pottypes6,'RowNames',label_subbasin_basin);
    
    tbl_percapita=table(subPot_GWh_percapita(:,1,ssel),subPot_GWh_percapita(:,2,ssel), subPot_GWh_percapita(:,3,ssel),...
        subPot_GWh_percapita(:,4,ssel), subPot_GWh_percapita(:,5,ssel), subPot_GWh_percapita(:,6,ssel),'VariableNames',pottypes6,'RowNames',label_subbasin_basin);    
    % 
    writetable(tbl_GWh,scfname,'Sheet',"GWh_yr",'WriteRowNames',true)
    writetable(tbl_perarea,scfname,'Sheet',"GWh_yr_area",'WriteRowNames',true) % per m^2
    writetable(tbl_percapita,scfname,'Sheet',"GWh_yr_capita",'WriteRowNames',true)
    end
end

%% Archive compiled runs
if write2mat
    theoryPot_ch=ch_basinEnergy;
    save(ofmatname, ...
        'pcsout','runnames', 'basindata', 'viz_compiled', 'theoryPot_ch',...
        'pcsouthaz', 'runnameshaz','pcsoutT2S', 'runnamesT2S','catchments_cl',...
        'subPot_GWh', 'subPot_num', 'subPot_GWh_perarea','subPot_GWh_percapita', ...
        'subPot_GWh_T2S')
    % For transport to R
    writetable(subPot_long,ofxls)
    disp('Data written to mat and xlsx')
end
% R.matlab is better/faster at handling long vectors than structures but
% saving to excel in matlab and loading that in R is faster
