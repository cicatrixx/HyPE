% Fixes order of saved theory pot data. plots basin totals for different river spacing.
% calculates and plots subbasin wise
% theory pot.

%Create final plot for: dQ vs dHP and Z vs HP
 % other plots are created in SubBasinWisePotential_Fut script
clear all
close all
addpath(genpath('Hydrus'))
run('myVarNames_Fut.m')
minQ=0.1; %in m3/s
miscCheck=0;
plotfigs=0;
% Load theoretical outputs
load(fullfile(rootof,"FutRuns_Figs_Analysis","fig_theory","FutTheoreticalPot_totals.mat"),...
    'fname','dorange','channel_GWh','Qannual_m3s_archive','channel_basinEnergy') % Q and HP data are both already re-ordered and masked

% Load subbasin boundaries
datapath=fullfile(pwd,"data","UI","data");
load(fullfile(datapath,'UI500m_ArcGIS.mat'),'catchments','basinlabels','outside','dem','channel')

disp("Loaded theory potential files")


%% Read correct order of models
modorder=readtable(fullfile(rootof,"FutRuns_Figs_Analysis","FutScen_Tracker.xlsx"),Sheet="GCMdetails", Range="A1:f13");
suffix=strcat("_LTavgs_",strrep(tframe,'-','_'));
myforder=[strcat(modorder.matlabname,suffix(1)); strcat(modorder.matlabname,suffix(2))];
newforder=zeros(1,length(myforder));

for f=1:length(myforder)
    newforder(f)=find(fname==myforder(f));
end
% keep historical at start
newforder=[1 newforder];
newfnames=fname(newforder');
channel_GWh_n=channel_GWh(:,newforder);
channel_basinEnergy_n=channel_basinEnergy(:,:,newforder);
Qannual_m3s_archive_n=Qannual_m3s_archive(:,:,newforder);
% prepare cleaned up names
ssp_gcm_names=extractBefore(newfnames(2:end),"_LTavgs_");
ssp_gcm_tf_names=["Historical"; strcat(ssp_gcm_names(1:12), "_",tframe(1)); strcat(ssp_gcm_names(13:24), "_",tframe(2))];
disp("Reorederd future scenarios")

%% GOOD: Plot hist vs fut total basin level theoretical potentail w bands for diff river segments lengths
% Setup start-end indices
rcpidx=[2:4:26]; % first index is historical
rcpselect=3*2:-1:1; % 3:-1:1;%
rcpsel=[1:3;4:6];

figure
for i =1:2
    subplot(1,2,i)
    hold all
    % Plot future simulations plot from reverse order so overlaps are
    % clearer
    for ssp=rcpsel(i,:)
        l(ssp)=plotEnvelopeMinMaxMean(channel_GWh_n(:,rcpidx(ssp):rcpidx(ssp+1)-1)/1000,c_ssp_main(ssp,:),0,dorange*cellsz_m/1000);
        scatter(dorange*cellsz_m/1000,channel_GWh_n(:,rcpidx(ssp):rcpidx(ssp+1)-1)/1000,5,c_ssp_main(ssp,:),'o');
    end
    l1=plot(dorange*cellsz_m/1000,channel_GWh_n(:,1)/1000,".-k","LineWidth",2); %,'DisplayName','Hindsight-based')%'Color',ccol(1,:),

    ylabel('Energy (TWh per year)','FontWeight','bold')
    xlabel('River segment length (km)','FontWeight','bold')
    applymyplotformat(sprintf("%s future: %s",tname(i), tframe(i)))
    ylim([800 2600])
end
l=legend([l1 l(rcpsel(i,:))],[ modelnames(2:end)],'Interpreter','none','Orientation','horizontal',Location='northoutside');

%% Evaluate basin and sub-basin wise theoretical potential
for f=1:length(fname)
    tmp_basinEnergy=channel_basinEnergy_n(:,:,f);
    for st=101:108 %basinlabels.basinIDs'
        selidxs=tmp_basinEnergy > 0 & catchments ==st;
        subPot_GWh(st-100,f) = sum(tmp_basinEnergy(selidxs),'all');
        subPot_num(st-100,:,f)=sum(selidxs,'all');
        subPot_GWh_per_m2(st-100,f) = subPot_GWh(st-100)/(sum(catchments ==st,'all')*cellsz_m^2);
    end
    % Add all basin total to row 9
    subPot_GWh(9,f) = sum(tmp_basinEnergy(tmp_basinEnergy > 0 & ~outside),'all');
    subPot_num(9,f) = sum(tmp_basinEnergy > 0 & ~outside,'all');

    subPot_GWh_per_m2(9,f) = subPot_GWh(st-100)/(sum(~outside,'all')*cellsz_m^2);
end
subPot_prctchange = (subPot_GWh-subPot_GWh(:,1))./subPot_GWh(:,1)*100;
subPot_names=[basinlabels.basinnames{1:8} "All Indus"];
mylim=[min(subPot_prctchange(1:8,:), [],'all') max(subPot_prctchange(1:8,:), [],'all')];

disp("Evaluated basin and sub-basinwise totals")

%% Write subbasin total and change in theoretical potential to file
oxlsfile=fullfile(rootoffut,'FutScen_Tracker.xlsx');
writetable(array2table(subPot_GWh,'VariableNames',ssp_gcm_tf_names','RowNames',subPot_names),oxlsfile,'Sheet','Theory_GWh','WriteRowNames', 1)
writetable(array2table(subPot_prctchange,'VariableNames',ssp_gcm_tf_names','RowNames',subPot_names),oxlsfile,'Sheet','Theory_prct','WriteRowNames', 1)
writetable(array2table(subPot_GWh_per_m2,'VariableNames',ssp_gcm_tf_names','RowNames',subPot_names),oxlsfile,'Sheet','Theory_GWh_m2','WriteRowNames', 1)
save(fullfile(rootoffut,'FutTheoreticalPot_subbasintotals.mat'), 'subPot_GWh','subPot_num','subPot_prctchange',"subPot_GWh_per_m2","subPot_names","ssp_gcm_tf_names")

disp("Theory pot totals written to excel and mat")

%% GOOD: Bar plot of change in basin-wide theory potential
figure
subplot(1,2,1)
bar(reshape(subPot_GWh(9,2:end)/1000,12,2)','EdgeAlpha',0)
yline(subPot_GWh(9,1)/1000,'color','red','LineWidth',1.5)
xline(1.5,'color','black','LineWidth',1.05) % tf splits

ylabel('Total energy (TWh/yr)','fontweight','bold')
applymyplotformat('Basin wide change in potential',c_ssp)
%text([1 2],2300*[1 1],strcat(tname,": ", tframe),'HorizontalAlignment','center') %,'fontweight','bold')
xticklabels(strcat(tname,": ", tframe))

subplot(1,2,2)
bar(reshape(subPot_prctchange(9,2:end),12,2)','EdgeAlpha',0)
xline(1.5,'color','black','LineWidth',1.05) % tf splits
ylabel("% Change in theoretical potential",'fontweight','bold')
applymyplotformat(' ',c_ssp)
%legend(ssp_gcm_names(1:12))
l=legend([repelem({' ',},24/2-4) cornernames],'NumColumns',3);
l.Title.String=strjoin(sspnames,'      ');
xticklabels(strcat(tname,": ", tframe))

%% Figs
if plotfigs
    %% GOOD: Bar plot of change in subbasin-wise potential
    %barplot is easier to read than scatter. so kept bar
    figure;
    subplot(3,4,1:4)
    bar(subPot_GWh(1:8,2:end)/1000,'EdgeAlpha',0)
    xline(1:8,'color','black','LineWidth',1.05) % tf splits
    ylabel('Total energy (TWh/yr)','fontweight','bold')
    xticklabels(basinlabels.basinnames)
    xtickangle(0)
    applymyplotformat('Sub-basin wise theoretical potential',c_ssp)
    l=legend([repelem({' ',},24/2-4) cornernames],'NumColumns',3);
    l.Title.String=strjoin(sspnames,'      ');

    for b=1:8
        subplot(3,4,b+4)
        bar(reshape(subPot_prctchange(b,2:end),12,[])','EdgeAlpha',0)
        %bar(reshape(subPot_prctchange(b,2:end),4,[])') this works but cant change color
        ylim([-25, 95])
        applymyplotformat(basinlabels.basinnames(b),c_ssp)
        %xline([4.5:4:24],'LineStyle',':') % ssp splits
        xline([1.5],'color','black','LineWidth',1.05) % tf splits
        xticklabels(strcat(tname,": ", tframe))
    end
    subplot(3,4,5)
    ylabel("% Change in theoretical potential",'fontweight','bold')
    %text([1 2],[91 91],strcat(tname,": ", tframe),'HorizontalAlignment','center')%,'fontweight','bold')

    %% GOOD: Spatial plot of subbasin wise change in theoretical potential
    catchments_cl=maskBasin(catchments,~outside);
    idata=subPot_prctchange; % basins are rowwise and scens are column wise
    subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005],[.08 .08],[.08 .01]); %[vgap hgap], marg_h -[low up],marg_w -[l r]

    figure; sgtitle(tf_full(1),'fontweight','bold');colormap(cbrewer2('Spectral') )
    c=1;sspi=1;
    for f=2:25
        if f==14; figure; sgtitle(tf_full(2),'fontweight','bold'); colormap(cbrewer2('Spectral') ); c=13; sspi=1; end
        spidx=f-c;
        subplot(3,4,spidx)
        tmpmap=changem(catchments_cl, subPot_prctchange(1:8,f), 101:108);
        h(f)=imagescnan(tmpmap);
        axis off
        caxis(mylim)
        % Annotations
        if spidx<=4;title(cornernames(spidx));end

        if ismember(spidx,[1:4:9])
            text(-500,1620/2,sspnames(sspi),'fontweight','bold','HorizontalAlignment','center','rotation',90);
            sspi=sspi+1;
        end
        if spidx==12;mycbar("Change in theoretical potential (%)",'southoutside'); end
    end

    %% Evaluate dchange in Q and HP at cell level
    nfuts=24;
    HP_GWh=maskBasin(channel_basinEnergy_n,~outside);
    HPchange_cellprct=(HP_GWh(:,:,2:end)-HP_GWh(:,:,1))./HP_GWh(:,:,1)*100;
    HPchange_cellprct_long=reshape(HPchange_cellprct,1620*2880,nfuts);
    myHPlim=[min(HPchange_cellprct,[],'all') max(HPchange_cellprct,[],'all')];

    Qannual_m3s=maskBasin(Qannual_m3s_archive_n,~outside);
    Qchange_cellprct=(Qannual_m3s(:,:,2:end)-Qannual_m3s(:,:,1))./Qannual_m3s(:,:,1)*100;
    Qchange_cellprct_long=reshape(Qchange_cellprct,1620*2880,nfuts);
    myQlim=[min(Qchange_cellprct,[],'all') max(Qchange_cellprct,[],'all')];
    disp("Evaluated cell level change in potential")

    %% Select HP, Q and dem for basin cells and plot linear regression between change in Q and HP only in the basin cells
    basinidx=find(outside==0);
    % Get data for just the basin cells
    tmpdem=dem(basinidx);
    for f=1:nfuts
        %get one future layer
        tmpQchange=Qchange_cellprct(:,:,f);
        tmpHPchange=HPchange_cellprct(:,:,f);
        %get basin cells for future layer
        tmpQ2_change(:,f)=tmpQchange(basinidx);
        tmpHP2_change(:,f)=tmpHPchange(basinidx);
        %diffQ_HP(f)=sum(tmpQ2(:,f)-tmpHP2(:,f),'omitnan');
    end
    %% Get HP energy at each elevation
    for f=1:nfuts+1 % incl historical
        tmpHP_GWh=HP_GWh(:,:,f);
        tmpHP2_GWh(:,f)=tmpHP_GWh(basinidx);
    end

     %% Plot elevation wise change in potential
    figure
    subplot(1,2,1)
    scatter(tmpHP2_change(:,1:12),tmpdem,20,'o','filled','MarkerFaceAlpha',0.5)
    applymyplotformat("Elevation vs dHP Mid",c_ssp)
    xlabel("Change in theoretical potential (%)")
    ylabel("Grid elevation (m)")

    subplot(1,2,2)
    scatter(tmpHP2_change(:,13:24),tmpdem,20,'o','filled','MarkerFaceAlpha',0.5)
    applymyplotformat("Elevation vs dHP Far",c_ssp)
    xlabel("Change in theoretical potential (%)")
    ylabel("Grid elevation (m)")
    l=  legend([repelem([" "],1,8), cornernames],'NumColumns',3);
    l.Title.String=strcat(rcpnames," ");

    %% FINAL: Evaluate and plot linear regression between change in Q and HP only in the basin cells
    fprintf("Difference between Q and HP basin cells is: %f \n",sum(sum(tmpQ2_change-tmpHP2_change,'omitnan')))
    % find indices for not nan cells in both Q and HP (these dont align across
    % all scenarios) there are fewer nans in Q cells
    valididx=find(~isnan(tmpQ2_change(:))& ~isnan(tmpHP2_change(:)));
    all(round(tmpQ2_change(valididx),3)==round(tmpHP2_change(valididx),3))

    % Linear fit between only notnan cells
    pp=polyfit(tmpQ2_change(valididx),tmpHP2_change(valididx),1);
    f=polyval(p,tmpQ2_change(valididx));

    figure;
    subplot(1,3,1)
    plot(tmpQ2_change(valididx),tmpHP2_change(valididx),'o',tmpQ2_change(valididx),f,'-')
    legend('All scenarios data','Linear fit x=y')
    xlabel("Change in annual Q (%)")
    ylabel("Change in theoretical potential (%)")
    applymyplotformat("Changes in each basin cell")
    axis square

    %% FINAL: Slow Plot elevation wise change in potential
    % manually setup format before plotting
    %figure
    %set(gcf, 'Position', get(0, 'Screensize'));
    subplot(1,3,2)
    for ff=[2:13,1,14:25,1]
        ll(ff)=plotCostCurve(tmpdem, tmpHP2_GWh(:,ff), ones(length(tmpdem),1),c_ssp26(ff,:),'-',myalpha,'',1);
        if ff==1 % shift to next subplot
            applymyplotformat(tf_full(1));
            ylabel("Elevation (m)")

            xlabel("Cumulative energy (TWh per year)")

            subplot(1,3,3)
            applymyplotformat(tf_full(2));
            ylabel("Elevation (m)")
            xlabel("Cumulative energy (TWh per year)")
        end
    end
    l=  legend(ll(1:13),histrcpcornernames);
    %l.Title.String=strcat(rcpnames," ");
   % savefig("tmp.fig")
end
%% Extras
if miscCheck
    %% Evaluate and plot dchange in basin theoretical potential over different river segments
    HP_prctchange = (channel_GWh_n-channel_GWh_n(:,1))./channel_GWh_n(:,1)*100;

    figure
    for i =1:2
        subplot(1,2,i)
        hold all
        % Plot future simulations plot from reverse order so overlaps are
        % clearer
        for ssp=rcpsel(i,:)
            l(ssp)=plotEnvelopeMinMaxMean(HP_prctchange(:,rcpidx(ssp):rcpidx(ssp+1)-1),c_ssp_main(ssp,:),0,dorange*cellsz_m/1000);
            scatter(dorange*cellsz_m/1000,HP_prctchange(:,rcpidx(ssp):rcpidx(ssp+1)-1),5,c_ssp_main(ssp,:),'o');

        end
        l1=plot(dorange*cellsz_m/1000,HP_prctchange(:,1),".-k","LineWidth",2); %,'DisplayName','Hindsight-based')%'Color',ccol(1,:),

        ylabel('Change in Theoretical Potential (%)','FontWeight','bold')
        xlabel('River segment length (km)','FontWeight','bold')
        applymyplotformat(sprintf("%s future: %s",tname(i), tframe(i)))
        %ylim([800 2600])
    end
    legend([l1 l(rcpsel(i,:))],["Historical", modelnames(3:end)],'Interpreter','none','Orientation','horizontal',Location='northoutside')

    %% SLOW: Plot Q and HP dchange in separate boxCHARTS
    figure
    subplot(2,1,1)
    boxchart(Qchange_cellprct_long,'MarkerStyle','none','Notch','on','JitterOutliers','on')  %'GroupByColor',fut_data.ssp,
    xline([4.5:4:24]) % ssp splits
    xline([12.5],'color','red','LineWidth',1.5) % tf splits
    applymyplotformat(strjoin(tframe,'                        '))
    xticks([])
    ylabel("Change in discharge (%)")
    ylim([-100 200])

    subplot(2,1,2)
    boxchart(HPchange_cellprct_long,MarkerStyle','none','JitterOutliers','on')  %'GroupByColor',fut_data.ssp,
    xline([4.5:4:24]) % ssp splits
    xline([12.5],'color','red','LineWidth',1.5) % tf splits
    set(gca,'XTick',1:24,'XTickLabel',fname(2:end),'ticklabelinterpreter','none','XTickLabelRotation',90)
    ylabel("Change in theoretical potential (%)")
    ylim([-100 200])
    applymyplotformat('')

    %% SLOW: Plot Q and HP dchange in the same boxPLOTS
    tbl1=table(repelem(1:nfuts,size(Qchange_cellprct_long,1))',Qchange_cellprct_long(:),repelem(["dQa"],numel(Qchange_cellprct_long))');
    tbl2=table(repelem(1:nfuts,size(HPchange_cellprct_long,1))',HPchange_cellprct_long(:),repelem(["dHP"],numel(HPchange_cellprct_long))');
    tbl=[tbl1;tbl2];
    tbl.Properties.VariableNames={'FutScen', 'Val', 'Vartype'};
    figure
    boxplot(tbl.Val,tbl.FutScen,'ColorGroup',tbl.Vartype,'whisker',inf, 'boxstyle','filled')%'PlotStyle','compact','OutlierSize',2,'Symbol','.')
    %,'label',{},'orientation','horizontal',ScensComp,'orientation','horizontal','FactorDirection','list','FullFactors','off','FactorSeparator',[1]); %'FactorGap',[10 ,.1],

    %% SLOW: Plot Q and HP dchange in the same boxCHARTS
    figure
    boxchart(tbl.Var1,tbl.Var2,'GroupByColor',tbl.Var3) %,'MarkerStyle','none')%'Notch','on',
    xline([4.5:4:24]) % ssp splits
    xline([12.5],'color','red','LineWidth',1.5) % tf splits
    set(gca,'XTick',1:24,'XTickLabel',fname(2:end),'ticklabelinterpreter','none','XTickLabelRotation',90)
    ylabel("% Change")
    legend(["Discharge" "Theoretical potential"])
    applymyplotformat('Changes in long term annual averages')
    ylim([-100 200])
    text([6 18],[180 180],strcat(tname,": ", tframe),'HorizontalAlignment','center','fontweight','bold')

    %% Plot change in cell level Q
    figure
    for i=2:13
        subplot(3,4,i-1)
        imagescnan(Qchange_cellprct(:,:,i))
        caxis(myQlim)
        %set(gca, 'ColorScale', 'log')
        %colorbar
    end
    cb = colorbar;
    %cb.Layout.Tile = 'east';
    m=1;
    for i=1:4:12
        subplot(3,4,i)
        ylabel(sspnames(m),'FontWeight','bold')
        m=m+1;
    end
    for i=1:4
        subplot(3,4,i)
        title(sprintf("Corner %d",i))
    end
    sgtitle("Change in Q (%)")

    %% Plot change in cell level H
    figure
    for i=2:13
        subplot(3,4,i-1)
        imagescnan(HPchange_cellprct(:,:,i))
        axis off
        caxis([-100 200])
        %set(gca, 'ColorScale', 'log')
        colorbar
    end
    cb = colorbar;
    %cb.Layout.Tile = 'east';
    m=1;
    for i=1:4:12
        subplot(3,4,i)
        ylabel(sspnames(m),'FontWeight','bold')
        m=m+1;
    end
    for i=1:4
        subplot(3,4,i)
        title(sprintf("Corner %d",i))
    end
    mycbar("Change in HP (%)")

    %% Plot subbasin-wise potential bars all in one
    %c_ssp=[c_ssp; brighten(c_ssp,0.4)];
    figure
    %total potential
    subplot(2,1,1)
    bar(subPot_GWh(:,2:end)/1000)
    %bar(diag(subPot_GWh/1000),'stacked') % for colorful bars
    ylabel('Total energy (TWh/yr)')
    xticklabels(basinlabels.basinnames)
    xtickangle(0)
    applymyplotformat(' ',c_ssp)

    % Change
    subplot(2,1,2)
    bar(subPot_prctchange(:,2:end))
    %bar(diag(subPot_GWh_perarea),'stacked') % for colorful bars
    ylabel('Change (%)')
    xticklabels(basinlabels.basinnames)
    xtickangle(0)
    applymyplotformat('',c_ssp)
    sgtitle("Sub-basin wise theoretical potential")

    l=legend([repelem({' ',},24/2-4) cornernames],'NumColumns',3);
    l.Title.String=strjoin(modelnames(3:end),'  ');

end
