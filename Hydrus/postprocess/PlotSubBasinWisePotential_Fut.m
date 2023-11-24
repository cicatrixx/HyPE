%% Get sub-basin wise totals for full and remaining scenarios

% Get sub-basin wise total potential, per area potential and per capita potential
% Cols are for theory, tech, fin, sust-tech, sust-fin, viz, exist. in some
% cases i removed the sust-tech entirely but not in all

% Also does comparison of future pop and demand vs potential

clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

ofmatname=fullfile(rootof,'MainFutScenarios.mat');
loadsaved=01;
idx=1;
plotallfigs=0;
write2mat=0;
write4cas=0;

%% Load catchment and Historic/fut population and demand data
datain=load(fullfile(rootf,sprintf('\\data\\%s\\Basin_UIB\\PantpeBasin_%d.mat',continent_in,nbasin)),'Regions', 'Pop','Countries');
cellsz_m=500;

load(fullfile(rootf,'data\UI\data\UI500m_ArcGIS.mat'), 'basinlabels','outside')
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

% Load future population and demand data
load(fullfile(rootf,'data','UI','data','futDemand_pop_energy.mat'),'futpopSubbasin','energySecReq_subbasin_TWh')

disp("Loaded basin info and pop+demand data")

%% Load or get subbasin wise data for all pot types
if loadsaved==0
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

    %% Get subbasin-wise THEORETICAL potential - all 25 scens
    theory=load(fullfile(rootoffut,'FutTheoreticalPot_subbasintotals.mat'), 'subPot_GWh', 'subPot_num');
    theory_subPot_GWh=theory.subPot_GWh;
    theory_subPot_num=theory.subPot_num;

    if plotallfigs
        %Plot theoretical and visualized potential
        figure;
        bar(theory_subPot_GWh,'FaceAlpha',0.5)
        hold all
        bar(viz_subPot_GWh,'FaceAlpha',0,'EdgeColor','r')
        bar(viz_subPot_GWh_exist,'FaceAlpha',0,'EdgeColor','g')
        legend(pottypes6([1,5,6]))
        xticklabels([basindata.basinnames{:},"All Indus"])
    end

    %% Get fnames of fut files in correct order of models (48 = 4cornersx 3SSP x 2TF x 2pottype)
    suffix=strrep(tframe,'-','_');
    suffix=[strcat("_",suffix(1),["_Tech_Fin" "_Sust_RiskAverse"]) strcat("_",suffix(2),["_Tech_Fin" "_Sust_RiskAverse"])]';
    myforder_mid=[];
    myforder_far=[];

    % Create neworder to have for each sspmodel in one place + mid first and
    % far after
    for f=1:height(modorder)
        myforder_mid=[myforder_mid ;strcat(modorder.matlabname(f),suffix(1:2))];
        myforder_far=[myforder_far ;strcat(modorder.matlabname(f),suffix(3:4))];
    end
    myforder=[myforder_mid;  myforder_far];

    % Get filepaths
    [fpaths, fnames]=path2fldrfiles(fullfile(rootf),strcat(r_prefix,"*"));
    % match order to new
    newforder=zeros(1,length(myforder));
    fnames_cl=strrep(strrep(strrep(fnames,"FutR103_",""),"_Full_Mixed",""),"_LTavgs","");
    for f=1:length(myforder)
        newforder(f)=find(fnames_cl==myforder(f));
    end
    newfnames=fnames_cl(newforder');
    newfpaths=fpaths(newforder);
    disp("Reorderd future scenarios")

    %% Load and compile historical potential portfolio
    hist=load(fullfile(rootof,'HistRuns_Figs_Analysis\Figs_trial\MainScenarios.mat'), 'runnames','pcsout');
    runnames=hist.runnames(5:6);
    pcsout{1}=hist.pcsout{5};
    pcsout{2}=hist.pcsout{6};
    disp("Loaded historical mixed portfolios")

    %% Loop through and compile future potential portfolio, remove heavy component data that i dont use so much
    runcnt=3;
    for fidx=1:length(newfpaths)
        % Load filename
        runnames{runcnt}=newfnames{fidx};
        pcsout{runcnt} = load(fullfile(newfpaths{fidx},sprintf('COEPOT_b%d_do.mat', nbasin)));
        pcsout{runcnt}.runname=newfnames{fidx};
        pcsout{runcnt}=rmfield(pcsout{runcnt}, {'RDDepth','RDCountry_id','RDRegion_id'} );
        runcnt=runcnt+1;
    end
    disp("Loaded future portfolios")

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

    disp("Assigned subbasin and country to each point")


    %% Calculate Sub-basin wise potential sums including existing
    potnames=["Theoretical", "Tech-Tech","Tech-Fin","Sust-Tech","Sust-Fin","Viz","Exist"];
    sc_idx=0;  % scenarios are layered in the 3Dmatrix
    for runid=1:numel(pcsout)
        fprintf("Compiling %s \n",runnames{runid})
        % Load and process techecon data
        mydata=pcsout{runid};
        if rem(runid,2)>0 % start new matrix sheet for odd runs and add theorypot to it
            sc_idx=sc_idx+1;
            % Add theory pot to start
            subPot_GWh(:,1,sc_idx)=theory_subPot_GWh(:,sc_idx);
            subPot_num(:,1,sc_idx)=theory_subPot_num(:,sc_idx);
            % Add vis pot to col 6
            subPot_GWh(:,6,sc_idx)=viz_subPot_GWh; % Full
            subPot_num(:,6,sc_idx)=viz_subPot_num;
            % Add exist pot to col 7
            subPot_GWh(:,7,sc_idx)=viz_subPot_GWh_exist; % Full
            subPot_num(:,7,sc_idx)=viz_subPot_num_exist;
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
                % Sust-Tech pot
                selidxs=subIDs ==st;
                subPot_GWh(st-100,4,sc_idx) = sum(GWhsel(selidxs));
                subPot_num(st-100,4,sc_idx)=sum(selidxs);
                % Sust-Fin pot
                selidxs=subIDs ==st & COEsel<=costlim;
                subPot_GWh(st-100,5,sc_idx) = sum(GWhsel(selidxs));
                subPot_num(st-100,5,sc_idx)=sum(selidxs);
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

            % Sust-Tech pot
            selidxs=subIDs>0;
            subPot_GWh(end,4,sc_idx) = sum(GWhsel(selidxs));
            subPot_num(end,4,sc_idx)=sum(selidxs);
            % Sust-Fin pot
            selidxs= COEsel<=costlim;
            subPot_GWh(end,5,sc_idx) = sum(GWhsel(selidxs));
            subPot_num(end,5,sc_idx)=sum(selidxs);
        end
    end

    disp("Calculated sub-basin wise totals")


    %% Get energy per area and per capita for subbasin
    for sc_idx=1:size(subPot_num,3)
        % For all subbasins
        for st=101:108
            subPot_GWh_perarea(st-100,:,sc_idx) =  subPot_GWh(st-100,:,sc_idx) /basindata.Area_m2(basindata.basinIDs==st);
            subPot_GWh_percapita_histpop(st-100,:,sc_idx) = subPot_GWh(st-100,:,sc_idx) /basindata.Pop_ppl(basindata.basinIDs==st);
        end
        % For all basin
        subPot_GWh_perarea(9,:,sc_idx) =  subPot_GWh(9,:,sc_idx)/basinArea_m2;
        subPot_GWh_percapita_histpop(9,:,sc_idx) = subPot_GWh(9,:,sc_idx)/basinPop_ppl; % this is only considering historical population
    end
    disp("Calculated per area and per capita vals")

    %% Eval future energy per capita and npop that can be sustained
    % Map population to potential years
    % Take avg pop for each time horizon. for mid: 2030-2050, for far: 2060-2080
    for selssp=1:3
        popSubbasin4terms(:,1,selssp) = mean(futpopSubbasin(1:9,3:5,selssp),2);
        popSubbasin4terms(:,2,selssp) = mean(futpopSubbasin(1:9,6:8,selssp),2);
    end

    % Eval % of fut pop that can be sustained and percapita energy
    subPot_PopSustained=subPot_GWh*1000/(energyreq_MWhpercapita); % convert GWH to MWH for division
    for selsssp=1:3
        %mid
        subPot_PopPrctofProjected(:,:,sspgrps{selsssp})=subPot_PopSustained(:,:,sspgrps{selsssp})./popSubbasin4terms(:,1,selssp)*100;
        %far
        subPot_PopPrctofProjected(:,:,sspgrps{selsssp+3})=subPot_PopSustained(:,:,sspgrps{selsssp+3})./popSubbasin4terms(:,2,selssp)*100 ;

        %mid
        subPot_MWh_percapita_futpop(:,:,sspgrps{selsssp})=subPot_GWh(:,:,sspgrps{selsssp})*1000./popSubbasin4terms(:,1,selssp);
        %far
        subPot_MWh_percapita_futpop(:,:,sspgrps{selsssp+3})=subPot_GWh(:,:,sspgrps{selsssp+3})*1000./popSubbasin4terms(:,2,selssp);
    end

    % For historical
    subPot_PopPrctofProjected(:,:,1)=subPot_PopSustained(:,:,1)./futpopSubbasin(1:9,1,1)*100;
    subPot_MWh_percapita_futpop(:,:,1)=subPot_GWh(:,:,1)*1000./futpopSubbasin(1:9,1,1);

    subPopPrctofProjectedChange=(subPot_PopSustained-subPot_PopSustained(:,:,1))./subPot_PopSustained(:,:,1)*100;
    subPopMWhPerCapitaPrctChange=(subPot_MWh_percapita_futpop-subPot_MWh_percapita_futpop(:,:,1))./subPot_MWh_percapita_futpop(:,:,1)*100;

    disp("Calculated future pop based params")

    %% Archive compiled runs
    % these subbasin numbers incl the sust-tech pot
    if write2mat

        save(ofmatname, 'newfnames', 'newfpaths',   ...
            'pcsout','runnames', 'basindata', 'viz_compiled',...
            'catchments_cl','subPot_GWh', 'subPot_num', 'subPot_GWh_perarea','subPot_GWh_percapita_histpop', 'potnames',...
            'subPot_PopPrctofProjected','subPopPrctofProjectedChange', 'subPot_MWh_percapita_futpop','subPopMWhPerCapitaPrctChange')

        disp('Data written to mat')
    end
else
    load(ofmatname, 'basindata','subPot_GWh','potnames','catchments_cl', 'subPot_PopPrctofProjected','subPopPrctofProjectedChange', 'subPot_MWh_percapita_futpop','subPopMWhPerCapitaPrctChange')
    disp("Loaded saved subPot data")
end
%% Calculate sub-basin wise change in potential
subPot_TWh=subPot_GWh/1000;
subPot_GWh_dchange=(subPot_GWh-subPot_GWh(:,:,1));
subPot_GWh_prct_change=(subPot_GWh-subPot_GWh(:,:,1))./subPot_GWh(:,:,1)*100;
basindata.basinnames(9)={'ALL INDUS'};
%skip the sust-tech pot in excel
selpot=[1 2 3 5 6 7 ];

%% Write to excel for CAS
if write4cas
    for ssel=1:size(subPot_GWh,3)
        tbl_GWh=array2table(subPot_GWh(:,selpot,ssel),'VariableNames',pottypes6,'RowNames',label_subbasin_basin);
        tbl_GWh_prctchange=array2table(subPot_GWh_prct_change(:,selpot,ssel),'VariableNames',pottypes6,'RowNames',label_subbasin_basin);
        tbl_MWh_percapita=array2table(subPot_MWh_percapita_futpop(:,selpot,ssel),'VariableNames',pottypes6,'RowNames',label_subbasin_basin);
        tbl_MWh_percapita_prctchange=array2table(subPopMWhPerCapitaPrctChange(:,selpot,ssel),'VariableNames',pottypes6,'RowNames',label_subbasin_basin);

        if ssel==1
            sbprefix="SubBasin_MediumTerm_";
        elseif ssel==14
            sbprefix="SubBasin_LongTerm_";
        end
        writetable(tbl_GWh,fullfile(rootoffut,strcat(sbprefix,'GWhperyr.xlsx')),'Sheet',strrep(histrcpcornernames{ssel},":","_"),'WriteRowNames',true)
        writetable(tbl_GWh_prctchange,fullfile(rootoffut,strcat(sbprefix,'GWhperyr_prctchange.xlsx')),'Sheet',strrep(histrcpcornernames{ssel},":","_"),'WriteRowNames',true)

        writetable(tbl_MWh_percapita,fullfile(rootoffut,strcat(sbprefix,'MWhpercapita.xlsx')),'Sheet',strrep(histrcpcornernames{ssel},":","_"),'WriteRowNames',true)
        writetable(tbl_MWh_percapita_prctchange,fullfile(rootoffut,strcat(sbprefix,'MWhpercapita_prctchange.xlsx')),'Sheet',strrep(histrcpcornernames{ssel},":","_"),'WriteRowNames',true)
    end
    disp('Data written to xlsx')
end

%% FINAL: Bar plot of total + spatial of far future for subbasin-wise potential - FOR ONE POTENTIAL ALL FUTS
%selscendata
for p=selpot(1:4)
    idataTWh=squeeze(subPot_TWh(:,p,:)); %2 and 4
    idataprct=squeeze(subPot_GWh_prct_change(:,p,:)); % basins are rowwise and scens are column wise - mid fut first, far fut second

    %barplot is easier to read than scatter. so kept bar
    figure;
    subplottight(5,4,1:4);
    b1=bar(idataTWh(1:8,2:end),'EdgeAlpha',0);
    hold on
    b2=bar(idataTWh(1:8,1),'FaceAlpha',0);
    xline(1:8,'color',mygraylines,'LineWidth',1.05) % tf splits
    xticklabels(basinlabels.basinnames)
    xtickangle(0)
    xlim([0.4,8.6])
    ylabel('Total energy (TWh/yr)','fontweight','bold')
    applymyplotformat(sprintf("Total %s potential",potnames{p}),c_ssp)
    l=legend([b1(13:24), b2],[repelem({' ',},24/2-4) cornernames 'Historical'],'NumColumns',4,'Location','northeast');
    l.Title.String=strjoin(rcpnames,'      ');
    % Add annotations
    addAlphaLabel(1,"outside")
    topy=ylim;
    text(1,max(idataTWh(1,:)),strcat("\leftarrow ", tname(1),"   " ,tname(2)," \rightarrow "),'HorizontalAlignment','center')%,'FontAngle','italic')
    
    % Shift downwards and increase height a bit
    h=gca;
    h.Position=[h.Position(1) h.Position(2)-.05 h.Position(3) h.Position(4)+.05];

    % Add spatial plots for far future
    seltf=2;
    mylim=[floor(min(idataprct(1:8,14:25),[],'all')) ceil(max(idataprct(1:8,14:25),[],'all'))];
    %mylim=[-15 95]; % fot theory pot
    tmpcmap=cbrewer2('RdBu',14); %skip the darkest and lightest shade
    colormap(tmpcmap(2:12,:))
    %colormap(viridis)
    c=9;sspi=1;
    for f=14:25
        subplottight(5,4,c)
        tmpmap=changem(catchments_cl, idataprct(1:8,f), 101:108);
        h(f)=imagescnan(tmpmap);
        axis off
        %colorbar
        caxis(mylim)
        % Add Annotations to top and left plots
        spidx=c-8;
        if spidx<=4
            title(cornernames(spidx));
        end
        if ismember(spidx,[1:4:9])
            text(-500,1620/2,rcpnames(sspi),'HorizontalAlignment','center','rotation',90,'fontweight','bold');
            sspi=sspi+1;
        end
        if spidx==9;
           cb= mycbar(sprintf("Change in %s potential in %s future (%%)",potnames{p},tf_full(seltf)),'southoutside');
            cb.Position=[0.09 0.06 0.9 0.015];    % shift it to bottom of fig
        end
        if spidx==1
                addAlphaLabel(2,"outside")
        end
        c=c+1;
    end
end

%% GOOD: Spatial plot of %change in  potential in mid and far future
for p=selpot(4)
    idata=squeeze(subPot_GWh_prct_change(:,p,:)); % basins are rowwise and scens are column wise - mid fut first, far fut second
    mylim=[floor(min(idata(1:8,:),[],'all')) ceil(max(idata(1:8,:),[],'all'))];
    figure;clf; sgtitle(tf_full(1),'fontweight','bold');
    colormap(viridis)
    c=1;sspi=1;
    for f=2:25
        if f==14; figure; sgtitle(tf_full(2),'fontweight','bold');     colormap(viridis); c=13; sspi=1; end
        spidx=f-c;
        subplottight(3,4,spidx)
        tmpmap=changem(catchments_cl, idata(1:8,f), 101:108);
        h(f)=imagescnan(tmpmap);
        axis off
        %colorbar
        caxis(mylim)
        % Annotations
        if spidx<=4;title(cornernames(spidx));end
        if ismember(spidx,[1:4:9])
            text(-500,1620/2,rcpnames(sspi),'fontweight','bold','HorizontalAlignment','center','rotation',90);
            sspi=sspi+1;
        end
        if spidx==9;mycbar(sprintf("Change in %s potential (%%)",potnames{p}),'southoutside');
        end
    end
end

%% GOOD: Bar plot of total and % change in subbasin-wise potential - FOR ONE POTENTIAL ALL FUTS
%barplot is easier to read than scatter. so kept bar
for p=selpot(4)
    idata=squeeze(subPot_GWh(:,p,:)); %2 and 4
    idata_prct_change=squeeze(subPot_GWh_prct_change(:,p,:)); %2 and 4

    figure;
    subplot(3,4,1:4)
    bar(idata(1:8,2:end)/1000,'EdgeAlpha',0)
    xline(1:8,'color','black','LineWidth',1.05) % tf splits
    ylabel('Total energy (TWh/yr)','fontweight','bold')
    xticklabels(basinlabels.basinnames)
    xtickangle(0)
    applymyplotformat(sprintf("Total %s potential",potnames{p}),c_ssp)
    l=legend([repelem({' ',},24/2-4) cornernames],'NumColumns',3);
    l.Title.String=strjoin(rcpnames,'      ');

    for b=1:8
        subplot(3,4,b+4)
        bar(reshape(idata_prct_change(b,2:end),12,[])','EdgeAlpha',0)
        %bar(reshape(subPot_prctchange(b,2:end),4,[])') this works but cant change color
        %ylim([-25, 95])
        applymyplotformat(basinlabels.basinnames(b),c_ssp)
        %xline([4.5:4:24],'LineStyle',':') % ssp splits
        xline([1.5],'color','black','LineWidth',1.05) % tf splits
        xticklabels(strcat(tname,": ", tframe))
    end
    subplot(3,4,5)
    ylabel("Change from historical (%%)",'fontweight','bold')
    %text([1 2],[91 91],strcat(tname,": ", tframe),'HorizontalAlignment','center')%,'fontweight','bold')
end

%% GOOD: Bar plot of only total subbasin-wise potential - FOR ONE POTENTIAL ALL FUTS
%barplot is easier to read than scatter. so kept bar
for p=selpot(1)
    idata=squeeze(subPot_GWh(:,p,:)); %2 and 4
    figure;
    b1=bar(idata(1:8,2:end)/1000,'EdgeAlpha',0);
    hold on
    b2=bar(idata(1:8,1)/1000,'FaceAlpha',0);
    xline(1:8,'color',mygraylines,'LineWidth',1.05) % tf splits
    xticklabels(basinlabels.basinnames)
    xtickangle(0)
    xlim([0.4,8.6])
    ylabel('Total energy (TWh/yr)','fontweight','bold')
    applymyplotformat(sprintf("Total %s potential",potnames{p}),c_ssp)
    l=legend([b1(13:24), b2],[repelem({' ',},24/2-4) cornernames 'Historical'],'NumColumns',4);
    l.Title.String=strjoin(rcpnames,'      ');
    topy=ylim;
    text(1,(topy(2)-100),strcat("\leftarrow ", tname(1),"   " ,tname(2)," \rightarrow "),'HorizontalAlignment','center')%,'FontAngle','italic')
end

%% GOOD: Bar plot only magnitude of change in subbasin-wise potential - FOR ONE POTENTIAL ALL FUTS
%barplot is easier to read than scatter. so kept bar
for p=selpot(4)
    idata=squeeze(subPot_GWh_dchange(:,p,:)); %2 and 4
    idata_prct_change=squeeze(subPot_GWh_prct_change(:,p,:)); %2 and 4

    figure;
    bar(idata(1:8,2:end)/1000,'EdgeAlpha',0)
    xline(1:8,'color',mygraylines,'LineWidth',1.05) % tf splits
    ylabel('\Delta TWh per yr','fontweight','bold')
    xticklabels(basinlabels.basinnames)
    xtickangle(0)
    applymyplotformat(sprintf("Total %s potential",potnames{p}),c_ssp)
    l=legend([repelem({' ',},24/2-4) cornernames],'NumColumns',3);
    l.Title.String=strjoin(rcpnames,'      ');
end

%% All-indus theory, tech, fin, sust potential, 4 corners separately - magnitude and change - only MF
selsubbas=9;
sel_sc=1:13; %only MF
sel_prct_change=squeeze(subPot_GWh_prct_change(selsubbas,selpot,sel_sc));
sel_subPot=squeeze(subPot_TWh(selsubbas,selpot,sel_sc));

figure
subplot(1,2,1)
plot(sel_subPot',"*")
set(gca,'xtick',sel_sc,'xticklabel',histrcpcornernames(sel_sc))
xline([1.5 5.5:4:13]) %,'LineWidth',1.5)
title("TWh per year")
%legend(pottypes3)
xlim([0.5 13.5])
grid on
subplot(1,2,2)
plot(sel_prct_change',"*")
set(gca,'xtick',sel_sc,'xticklabel',histrcpcornernames(sel_sc))
xline([1.5 5.5:4:13]) %,'LineWidth',1.5)
yline(0,'Color','b')
title("% Change from historical")
sgtitle(strcat(basindata.basinnames(selsubbas )," for Near future"))
%legend(pottypes3)
xlim([1.5 13.5])
legend(pottypes6([1:4,5,6]),'Location','southoutside','NumColumns',2)
grid on

%% FINAL: Sub-basin % change in potentials 4 corners separately scatter w | - both MF and FF
% + or . symbols dont work so well
figure
for selsubbas=1:9
    nexttile
    plot([1 repelem([2:7],4)]-.35,squeeze(subPot_GWh_prct_change(selsubbas,selpot(1),:))',"|")
    hold all
    plot([1 repelem([2:7],4)]-.15,squeeze(subPot_GWh_prct_change(selsubbas,selpot(2),:))',"|")
    plot([1 repelem([2:7],4)]+.15,squeeze(subPot_GWh_prct_change(selsubbas,selpot(3),:))',"|")
    plot([1 repelem([2:7],4)]+.35,squeeze(subPot_GWh_prct_change(selsubbas,selpot(4),:))',"|")

    set(gca,'xtick',1:7,'xticklabel',["Historical" tfrcpnames]) %Xticklength
    xline([1.5 4.5]) %,'LineWidth',1.5)
    title(basindata.basinnames(selsubbas ))
    xlim([1.5 7.5])
    if selsubbas<9
        ylim([-50 150])
    end
    grid on
    %set(gca,'YScale', 'log')
    %colororder(hsv(4))
    %colororder(cmap8_wong([1 4 5 3 6 7],:))
end
set(gca,'YColor', [165, 42, 42]/255); %brown

legend(pottypes6([1:4]),'Location','southoutside','NumColumns',2)
sgtitle("% change in Subbasin Potentials")

%% GOOD: Sub-basin total potentials 4 corners separately scatter w | - both MF and FF
% + or . symbols dont work so well
figure
for selsubbas=1:9
    nexttile
    hold all

    % Shift pots left and right
    p(1)=plot([1 repelem([2:7],4)]-.35,squeeze(subPot_TWh(selsubbas,selpot(1),:))',"|", DisplayName=pottypes6{1});
    p(2)=plot([1 repelem([2:7],4)]-.15,squeeze(subPot_TWh(selsubbas,selpot(2),:))',"|", DisplayName=pottypes6{2});
    p(3)=plot([1 repelem([2:7],4)]+.15,squeeze(subPot_TWh(selsubbas,selpot(3),:))',"|", DisplayName=pottypes6{3});
    p(4)=plot([1 repelem([2:7],4)]+.35,squeeze(subPot_TWh(selsubbas,selpot(4),:))',"|", DisplayName=pottypes6{4});

    set(gca,'xtick',1:7,'xticklabel',["Historical" tfrcpnames]) %Xticklength
    xline([1.5 4.5]) %,'LineWidth',1.5)
    applymyplotformat(basindata.basinnames(selsubbas ))
    xlim([1.5 7.5])
    if selsubbas<9
        ylim([0 100])
        %    ylim([0.1 1100]) % if theory pot incl

    end
    grid on
    %set(gca,'YScale', 'log')
    %colororder(hsv(4))
    %colororder(cmap8_wong([1 4 5 3 6 7],:))
end
set(gca,'YColor', [165, 42, 42]/255); %brown

legend(p,'Location','southoutside','NumColumns',2)
sgtitle("TWh per year")

%% Sub-basin total all 6 pots, 4 corners separately - only magnitude %only MF
sel_sc=1:13;
figure
for selsubbas=1:9
    nexttile
    sel_subPot=squeeze(subPot_TWh(selsubbas,selpot,sel_sc));
    plot(sel_subPot',"*")
    set(gca,'xtick',sel_sc,'xticklabel',histrcpcornernames(sel_sc))
    xline([1.5 5.5:4:13]) %,'LineWidth',1.5)
    ylabel("TWh per year")
    xlim([0.5 13.5])
    grid on
    title(basindata.basinnames(selsubbas))
end
legend(pottypes6([1:4,5,6]),'Location','southoutside','NumColumns',2)
grid on
sgtitle("for Mid future")

%% Sub-basin sample, all 6 pots, all corners in one line - only MF
sel_sc=1:13; %only MF
figure
for selsubbas=1:9
    nexttile
    %subplot(1,2,1)
    plot([1 repelem([2 3 4],4)],squeeze(subPot_TWh(selsubbas,selpot,sel_sc))',"*")
    set(gca,'xtick',1:4,'xticklabel',["Hist" rcpnames])
    %xline([1.5 5.5:4:13]) %,'LineWidth',1.5)
    title(basindata.basinnames(selsubbas ))
    xlim([0.5 4.5])
    ylim([10^(-1) 2000])
    grid on
    %set(gca,'YScale', 'log')
end
legend(pottypes6([1:4,5,6]),'Location','southoutside','NumColumns',2)
sgtitle("TWh/year for Mid fut")

%% GOOD: Sub-basin potentials all corners in one line - both MF and FF
figure
for selsubbas=1:9
    nexttile
    %subplot(1,2,1)
    plot([1 repelem([2:7],4)],squeeze(subPot_TWh(selsubbas,selpot(2:4),:))',"*")
    hold all
    plot([1 repelem([2:7],4)],squeeze(subPot_TWh(selsubbas,selpot(6),:))',"x")

    set(gca,'xtick',1:7,'xticklabel',["Historical" rcpnames],'XTickLabelRotation',90)
    xline([1.5 4.5]) %,'LineWidth',1.5)
    title(basindata.basinnames(selsubbas ))
    xlim([0.5 7.5])
    if selsubbas<9
        ylim([0 100])
        %ylim([0.1 1100]) % if theory pot incl

    end
    grid on
    %set(gca,'YScale', 'log')
    colororder(cmap8_wong_gray([1 5 3 7],:))
end
set(gca,'YColor', [165, 42, 42]/255); %brown

legend(pottypes6([2:4,6]),'Location','southoutside','NumColumns',2)
sgtitle("Subbasin Potentials in TWh/year")


%% GOOD: Sub-basin dchange in potentials all corners in one line - both MF and FF
figure
for selsubbas=1:9
    nexttile
    plot([1 repelem([2:7],4)],squeeze(subPot_GWh_dchange(selsubbas,selpot(1),:))',"o")
    hold all
    plot([1 repelem([2:7],4)],squeeze(subPot_GWh_dchange(selsubbas,selpot(2),:))',"*")
    plot([1 repelem([2:7],4)],squeeze(subPot_GWh_dchange(selsubbas,selpot(3),:))',"+")
    plot([1 repelem([2:7],4)],squeeze(subPot_GWh_dchange(selsubbas,selpot(4),:))',".")

    %plot([1 repelem([2:7],4)],squeeze(subPot_GWh_dchange(selsubbas,selpot(5:6),:))',"x")

    set(gca,'xtick',1:7,'xticklabel',["Historical" tfrcpnames],'XTickLabelRotation',90)
    xline([1.5 4.5]) %,'LineWidth',1.5)
    title(basindata.basinnames(selsubbas ))
    if selsubbas<9
        %ylim([0.1 1100])
        % xlim([0.5 7.5])

    end
    grid on
    %set(gca,'YScale', 'log')
    colororder(cmap8_wong_gray([1 4 5 3 6 7],:))
end
set(gca,'YColor', [165, 42, 42]/255); %brown

legend(pottypes6([1:4]),'Location','southoutside','NumColumns',2)
sgtitle("TWh change in Subbasin Potentials")

%% Calculate min/max mean for each scenario for each subbasin
% For total pot
for selsubbas=1:9
    pp=1; % because selpot is not serial
    for p=selpot(1:4)
        idata=squeeze(subPot_TWh(selsubbas,p,:)); % basins are rowwise and scens are column wise - mid fut first, far fut second
        rcpgroups=[0 repelem([1:6],4)]';

        % Save mean in 3d mat with rcp x pottype x subbasin
        idata_byrcp(:,:,pp)=reshape(idata(2:end),[4  6]);
        idata_mean(:,pp,selsubbas)=groupsummary(idata,rcpgroups,'mean') ;
        idata_neg(:,pp,selsubbas)= groupsummary(idata,rcpgroups,'min') - idata_mean(:,pp,selsubbas);
        idata_pos(:,pp,selsubbas)= groupsummary(idata,rcpgroups,'max')- idata_mean(:,pp,selsubbas);
        pp=pp+1;
    end
    subbasdata_byrcp{selsubbas}=idata_byrcp;
end

% For % change
for selsubbas=1:9
    pp=1; % because selpot is not serial
    for p=selpot(1:4)
        idata=squeeze(subPot_GWh_prct_change(selsubbas,p,:)); % basins are rowwise and scens are column wise - mid fut first, far fut second
        rcpgroups=[0 repelem([1:6],4)]';

        % Save mean in 3d mat with rcp x pottype x subbasin
        idata_byrcp_prct(:,:,pp)=[repelem(idata(1),4)' reshape(idata(2:end),[4  6])];
        idata_mean_prct(:,pp,selsubbas)=groupsummary(idata,rcpgroups,'mean') ;
        idata_neg_prct(:,pp,selsubbas)=groupsummary(idata,rcpgroups,'min') -idata_mean_prct(:,pp,selsubbas);
        idata_pos_prct(:,pp,selsubbas)=groupsummary(idata,rcpgroups,'max')- idata_mean_prct(:,pp,selsubbas);
        pp=pp+1;
    end
    subbasdata_byrcp_prct{selsubbas}=idata_byrcp_prct;
end

nscens=size(idata_mean,1);
npots=size(idata_mean(:,:,selsubbas),2);
disp("Reshaped data for rcp vs corner")


%% GOOD: Range plot for sub-basin total potentials - incl theoretical
c_sel= c_pot5; %hsv(4); %cmap8_wong_gray;

xcenter=1:nscens;
xshifted=[ xcenter-.35
    xcenter-.15
    xcenter+.15
    xcenter+.35];
figure
for selsubbas=1:9
    nexttile
    hold all
    for pp=1:npots
        ee(pp)=errorbar(xshifted(pp,:),idata_mean(:,pp,selsubbas), idata_neg(:,pp,selsubbas), idata_pos(:,pp,selsubbas),'_', ...
            'MarkerSize',0.1,'CapSize',0,'LineWidth',4.5,'Color',c_sel(pp,:),'DisplayName',pottypes6{pp});
        plot(xshifted(pp,:),idata_mean(:,pp,selsubbas),".",'Color',mygraylines)
    end
    set(gca,'xtick',1:7,'xticklabel',["Historical" tfrcpnames]) %Xticklength
    xline([1.5 4.5]) %,'LineWidth',1.5)
    applymyplotformat(basindata.basinnames(selsubbas ))
    %xlim([1.5 7.5])
    if selsubbas<9
        %ylim([0 100])
        ylim([0.1 1100]) % if theory pot incl
    end
    grid on
end
set(gca,'YColor', [165, 42, 42]/255); %brown

legend(ee,'Location','southoutside','NumColumns',2)
sgtitle("TWh per year")


%% FINAL: Range plot for sub-basin total potentials - NO theoretical
c_sel= c_pot5; %hsv(4); %cmap8_wong_gray;
xcenter=1:7;
xshifted=[ xcenter-.35
    xcenter-.15
    xcenter+.15
    xcenter+.35];
figure
for selsubbas=[9 1:8]
    nexttile
    hold all

    for pp=2:npots
        ee1(pp-1)=errorbar(xshifted(pp,:),idata_mean(:,pp,selsubbas), idata_neg(:,pp,selsubbas), idata_pos(:,pp,selsubbas),'_', ...
            'MarkerSize',0.1,'CapSize',0,'LineWidth',4.5,'Color',c_sel(pp,:),'DisplayName',pottypes6{pp});
        ee2=plot(xshifted(pp,:),idata_mean(:,pp,selsubbas),"_k",'DisplayName','Mean of corners'); %,'Color',mygraylines)
    end
    set(gca,'xtick',1:7,'xticklabel',["Historical" tfrcpnames]) %Xticklength
    xline([1.5 4.5]); %,'LineWidth',1.5)
    % Add existing pot as a line
    ee3=yline(subPot_TWh(selsubbas,7,1),'Linestyle',':','Color','red','LineWidth',1.25   ,'DisplayName','Existing');
    applymyplotformat(basindata.basinnames(selsubbas ))
    set(gca,'XGrid','off','YMinorGrid','on' )
    xlim([0.5 7.5])
    ylabel("TWh per year")

    if selsubbas<9
        ylim([0 100])
        % ylim([0.1 1100]) % if theory pot incl
    end

    if selsubbas==9
        set(gca,'YColor',extraaxis ,'XColor',extraaxis );
    end

end
legend([ee1 ee3 ee2(1)],'Location','southoutside','NumColumns', 3, 'Orientation', 'horizontal')
sgtitle("TWh per year")

%% FINAL: Range plot for sub-basin % change in potentials - WITH theoretical
c_sel= c_pot5; %hsv(4); %cmap8_wong_gray;
xcenter=1:7;
xshifted=[ xcenter-.35
    xcenter-.15
    xcenter+.15
    xcenter+.35];
figure
for selsubbas=[9 1:8]
    nexttile
    hold all

    for pp=1:npots
        ee1(pp)=errorbar(xshifted(pp,:),idata_mean_prct(:,pp,selsubbas), idata_neg_prct(:,pp,selsubbas), idata_pos_prct(:,pp,selsubbas),'_', ...
            'MarkerSize',0.1,'CapSize',0,'LineWidth',4,'Color',c_sel(pp,:),'DisplayName',pottypes6{pp});
        ee2=plot(xshifted(pp,:),idata_mean_prct(:,pp,selsubbas),"_k",'DisplayName','Mean of corners','MarkerSize',4 ); %,'Color',mygraylines)
    end
    set(gca,'xtick',1:7,'xticklabel',["Historical" tfrcpnames]) %Xticklength
    xline([1.5 4.5]); %,'LineWidth',1.5)
    applymyplotformat(label_subbasin_basin(selsubbas ))
    xlim([1.5 7.5]) % can skip the historical one
    ylabel("% change")
    set(gca,'XGrid','off','YMinorGrid','on' )

    if selsubbas<9
        ylim([-50 150])
        yticks(-50:25:150)
    end
    if selsubbas==9
        set(gca,'YColor',extraaxis ,'XColor',extraaxis );
    end
end

legend([ee1 ee2(1)],'Location','southoutside','NumColumns', 3, 'Orientation', 'horizontal')
sgtitle("% change in potential")

%% >>>>>>>>>GOOD: Spatial plot of per capita energy availability for SUST FAR
for p=selpot(4)
    idata=squeeze(subPot_MWh_percapita_futpop(1:8,p,:)); % basins  are rowwise and scens are column wise - mid fut first, far fut second
    mylim=[floor(min(idata,[],'all')) ceil(max(idata,[],'all'))];
    figure;clf; sgtitle(tf_full(1),'fontweight','bold');colormap(cbrewer2('Spectral') )
    c=1;sspi=1;
    for f=2:25
        if f==14; figure; sgtitle(tf_full(2),'fontweight','bold'); colormap(cbrewer2('Spectral') ); c=13; sspi=1; end
        spidx=f-c;
        subplottight(3,4,spidx)
        tmpmap=changem(catchments_cl, idata(1:8,f), 101:108);
        h(f)=imagescnan(tmpmap);
        axis off
        % colorbar
        %colormap(parula(16))
        caxis(mylim)

        %set(gca,'ColorScale','log')
        % Annotations
        if spidx<=4;title(cornernames(spidx));end

        if ismember(spidx,[1:4:9])
            text(-500,1620/2,sspnames(sspi),'fontweight','bold','HorizontalAlignment','center','rotation',90);
            sspi=sspi+1;
        end
        if spidx==9
            cb=mycbar(sprintf("Future MWh per capita per year for %s potential",potnames{p}),'southoutside');
            %set(cb,'ytick',[0:50:3000])
        end
    end
    cc=colormap;
    figure
    tmpmap=changem(catchments_cl, idata(1:8,1), 101:108);
    imagescnan(tmpmap);
    axis off
    caxis(mylim)
    colormap(cc)
    cb=mycbar(sprintf("Historical MWh per capita per year for %s potential",potnames{p}),'southoutside');
end

%% GOOD: Spatial plot of % change in per capita energy availability
for p=selpot(4)
    idata=squeeze(subPopMWhPerCapitaPrctChange(1:8,p,:)); % basins  are rowwise and scens are column wise - mid fut first, far fut second
    mylim=[floor(min(idata,[],'all')) ceil(max(idata,[],'all'))];
    figure;clf; sgtitle(tf_full(1),'fontweight','bold');colormap(cbrewer2('Spectral') )
    c=1;sspi=1;
    for f=2:25
        if f==14; figure; sgtitle(tf_full(2),'fontweight','bold'); colormap(cbrewer2('Spectral') ); c=13; sspi=1; end
        spidx=f-c;
        subplottight(4,4,spidx)
        tmpmap=changem(catchments_cl, idata(1:8,f), 101:108);
        h(f)=imagescnan(tmpmap);
        axis off
        % colorbar
        %colormap(parula(16))
        caxis(mylim)

        %set(gca,'ColorScale','log')
        % Annotations
        if spidx<=4;title(cornernames(spidx));end

        if ismember(spidx,[1:4:9])
            text(-500,1620/2,sspnames(sspi),'fontweight','bold','HorizontalAlignment','center','rotation',90);
            sspi=sspi+1;
        end
        if spidx==9
            cb=mycbar(sprintf("%% change in Future MWh per capita per year for %s potential",potnames{p}),'southoutside');
            %set(cb,'ytick',[0:50:3000])
        end
    end

end

%% FINAL: Bar plot of change in subbasin-wise per capita potential
%barplot is easier to read than scatter. so kept bar
clear subplot
for p=selpot([2,4])
    idata=squeeze(subPot_MWh_percapita_futpop(:,p,:)); %2 and 4
    idata_prct_change=squeeze(subPopMWhPerCapitaPrctChange(:,p,:)); %2 and 4

    figure;
    subplot(3,4,1:4)
    b1=bar(idata(1:8,2:end),'EdgeAlpha',0);
    hold on
    b2=bar(idata(1:8,1),'FaceAlpha',0);
    xline(1:8,'color',mygraylines,'LineWidth',1.05) % tf splits
    xticklabels(basinlabels.basinnames)
    xtickangle(0)
    ylabel('MWh per capita per year','fontweight','bold')
    applymyplotformat(sprintf("Energy availability under %s ",potnames{p}),c_ssp)
    % Add min energy threshold
    b3=yline(energyreq_MWhpercapita,':r');

    l=legend([b1(13:24), b2, b3],[repelem({' ',},24/2-4) cornernames 'Historical' 'Min electricity req'],'NumColumns',4);
    l.Title.String=strjoin(rcpnames,'      ');

    for b=1:8
        subplot(3,4,b+4)
        bar(reshape(idata_prct_change(b,2:end),12,[])','EdgeAlpha',0)
        applymyplotformat(basindata.basinnames(b),c_ssp)
        %xline([4.5:4:24],'LineStyle',':') % ssp splits
        xline([1.5],'color','black','LineWidth',1.05) % tf splits
        xticklabels(strcat(tname,": ", tframe))
    end
    subplot(3,4,5)
    ylabel("% Change from historical",'fontweight','bold')
    %text([1 2],[91 91],strcat(tname,": ", tframe),'HorizontalAlignment','center')%,'fontweight','bold')
end


%% GOOD: Spatial plot of change in % of projected pop that can be supported
for p=selpot(4)
    idata=squeeze(subPot_PopPrctofProjected(1:8,p,:)); % basins  are rowwise and scens are column wise - mid fut first, far fut second
    mylim=[floor(min(idata,[],'all')) ceil(max(idata,[],'all'))];
    figure;clf; sgtitle(tf_full(1),'fontweight','bold');colormap(cbrewer2('Spectral') )
    c=1;sspi=1;
    for f=2:25
        if f==14; figure; sgtitle(tf_full(2),'fontweight','bold'); colormap(cbrewer2('Spectral') ); c=13; sspi=1; end
        spidx=f-c;
        subplottight(3,4,spidx)
        tmpmap=changem(catchments_cl, idata(1:8,f), 101:108);
        h(f)=imagescnan(tmpmap);
        axis off
        % colorbar
        %colormap(parula(16))
        caxis(mylim)

        %set(gca,'ColorScale','log')
        % Annotations
        if spidx<=4;title(cornernames(spidx));end

        if ismember(spidx,[1:4:9])
            text(-500,1620/2,sspnames(sspi),'fontweight','bold','HorizontalAlignment','center','rotation',90);
            sspi=sspi+1;
        end
        if spidx==9
            cb=mycbar(sprintf("Percentage of sub-basin population that can be\n supplied %.2f MWh electricity annually by %s potential",energyreq_MWhpercapita,potnames{p}),'southoutside');
            %set(cb,'ytick',[0:50:3000])
        end
    end
    cc=colormap;
    figure
    tmpmap=changem(catchments_cl, idata(1:8,1), 101:108);
    imagescnan(tmpmap);
    caxis(mylim)
    colormap(cc)
end

%% GOOD: Spatial plot of change in % of projected pop that can be supported CATEGORICAL
mycats=[0 100 200 400 800 1000 Inf];
ncats=length(mycats)-1;
mycatlabels=[compose("<=%d%%",mycats(2:ncats)), ">1000"];
cmap = lines(ncats);

for p=selpot(4)
    idata0=squeeze(subPot_PopPrctofProjected(1:8,p,:)); % basins  are rowwise and scens are column wise - mid fut first, far fut second
    idata=discretize(idata0,mycats);
    mylim=[floor(min(idata,[],'all')) ceil(max(idata,[],'all'))];
    figure;clf; sgtitle(tf_full(1),'fontweight','bold');colormap(cmap);

    c=1;sspi=1;
    for f=2:25
        if f==14; figure; sgtitle(tf_full(2),'fontweight','bold');colormap(cmap); c=13; sspi=1; end
        spidx=f-c;
        subplottight(3,4,spidx)
        tmpmap=changem(catchments_cl, idata(1:8,f), 101:108);
        h(f)=imagescnan(tmpmap);
        axis off
        % colorbar
        %colormap(parula(16))
        caxis(mylim)

        %set(gca,'ColorScale','log')
        % Annotations
        if spidx<=4;title(cornernames(spidx));end

        if ismember(spidx,[1:4:9])
            text(-500,1620/2,sspnames(sspi),'fontweight','bold','HorizontalAlignment','center','rotation',90);
            sspi=sspi+1;
        end
        if spidx==9
            cb=mycbar(sprintf("Percentage of sub-basin population that can be\n supplied %.2f MWh electricity annually by %s potential",energyreq_MWhpercapita,potnames{p}),'southoutside');
            % Create a colorbar and set its properties
            cb.Ticks = (1:ncats);
            cb.TickLabels = mycatlabels;
            %set(cb,'ytick',[0:50:3000])
        end
    end

    figure
    tmpmap=changem(catchments_cl, idata(1:8,1), 101:108);
    imagescnan(tmpmap);
    caxis(mylim)
    colormap(cmap)
    title("historical")
    colormap(cmap)
    cb=mycbar(sprintf("Percentage of sub-basin population that can be\n supplied %.2f MWh electricity annually by %s potential",energyreq_MWhpercapita,potnames{p}),'southoutside');
    % Create a colorbar and set its properties
    cb.Ticks = 1:ncats;
    cb.TickLabels = mycatlabels;
    %set(cb,'ytick',[0:50:3000])
end

%% GOOD: Bar plot of change in subbasin-wise population supported
%barplot is easier to read than scatter. so kept bar
for p=selpot([2,4])
    idata=squeeze(subPot_PopPrctofProjected(:,p,:)); %2 and 4
    idata_prct_change=squeeze(subPopPrctofProjectedChange(:,p,:)); %2 and 4

    figure;
    subplot(3,4,1:4)
    b1=bar(idata(1:8,2:end),'EdgeAlpha',0);
    hold on
    b2=bar(idata(1:8,1),'FaceAlpha',0);
    xline(1:8,'color',mygraylines,'LineWidth',1.05) % tf splits
    xticklabels(basinlabels.basinnames)
    xtickangle(0)
    ylabel('% of sub-basin pop','fontweight','bold')
    applymyplotformat(sprintf("Total %s of sub-basin popln",potnames{p}),c_ssp)

    l=legend([repelem({' ',},24/2-4) cornernames],'NumColumns',3,'Location','bestoutside');
    l.Title.String=strjoin(rcpnames,'      ');

    for b=1:8
        subplot(3,4,b+4)
        bar(reshape(idata_prct_change(b,2:end),12,[])','EdgeAlpha',0)
        %bar(reshape(subPot_prctchange(b,2:end),4,[])') this works but cant change color
        %ylim([-25, 95])
        applymyplotformat(basindata.basinnames(b),c_ssp)
        %xline([4.5:4:24],'LineStyle',':') % ssp splits
        xline([1.5],'color','black','LineWidth',1.05) % tf splits
        xticklabels(strcat(tname,": ", tframe))
    end
    subplot(3,4,5)
    ylabel("% Change from historical",'fontweight','bold')
    %text([1 2],[91 91],strcat(tname,": ", tframe),'HorizontalAlignment','center')%,'fontweight','bold')
end




