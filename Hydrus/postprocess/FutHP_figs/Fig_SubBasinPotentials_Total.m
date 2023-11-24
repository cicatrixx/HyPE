% Fig_SubBasinPotentials_Total
% Created By    : Sanita Dhaubanjar
% Created For	: SustainIndus project WP2 Paper publication
%=========================
% Code for loading pre-processed output files from Hydrus model to create
% final plots used for publication
% Plot sub-basin wise bar for total potential and spatial map for % change in potential for far
% future

clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

%% Load saved subbasin data for fut potentials
ofmatname=fullfile(rootof,"FutFigs_cl0","FutFigs_cl1","FigData",'MainFutScenarios.mat');
load(ofmatname,'subPot_GWh','potnames','catchments_cl')
disp("Loaded saved subPot data")
potnames_cl5=[potnames(1) newpots4cl];

%% Calculate sub-basin wise change in potential
subPot_TWh=subPot_GWh/1000;
subPot_prct_change=(subPot_GWh-subPot_GWh(:,:,1))./subPot_GWh(:,:,1)*100;
%skip the sust-tech pot in excel
selpot=[1 2 3 5 6 7 ];

%% FINAL: Bar plot of total + spatial of far future for subbasin-wise potential - FOR ONE POTENTIAL ALL FUTS
%selscendata
for p=selpot(4) % 4=sustainable potential
    idataTWh=squeeze(subPot_TWh(:,p,:)); %2 and 4
    idataprct=squeeze(subPot_prct_change(:,p,:)); % basins are rowwise and scens are column wise - mid fut first, far fut second

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
    ylabel(energylabel,'fontweight','bold','FontSize',11)
    applymyplotformat(sprintf("Total %s potential",potnames_cl5{p}),c_ssp)
    l=legend([b1(13:24), b2],[repelem({' ',},24/2-4) cornernames 'Historical'],'NumColumns',4,'Location','northeast');
    l.Title.String=strjoin(rcpnames,'      ');
    % Add annotations
    addAlphaLabel(1,"outside")
    topy=ylim;
    %text(1,max(idataTWh(1,:)),strcat("\leftarrow ", tname(1),"   " ,tname(2)," \rightarrow "),'HorizontalAlignment','center')%,'FontAngle','italic')
    text(1,floor(topy(2)),strcat(tname(1),"      " ,tname(2)),'HorizontalAlignment','center')%,'FontAngle','italic')

    % Shift downwards and increase height a bit
    h=gca;
    h.Position=[h.Position(1) h.Position(2)-.05 h.Position(3) h.Position(4)+.05];

    % Add spatial plots for far future
    seltf=2;
    %get min-max for far future scenarios to set limits for all cmaps
    mylim=[floor(min(idataprct(1:8,14:25),[],'all')) ceil(max(idataprct(1:8,14:25),[],'all'))];
    if p==selpot(4) % For sustainable pot
        mylim=[-30 150];
        %colormap(cbrewer2('RdYlBu'))
        tmpcmap=[brighten(flipud(cbrewer2('YlOrRd',50)),0.1);  cbrewer2('YlGnBu',185)];
        colormap(tmpcmap([10:45,71:end],:))
        %colormap(xx)
    end
    c=9;sspi=1;
    for f=14:25
        subplottight(5,4,c)
        tmpmap=changem(catchments_cl, idataprct(1:8,f), 101:108);
        h(f)=imagescnan(tmpmap);
        axis off
        caxis(mylim)

        % add Annotations to top and left plots
        spidx=c-8;
        if spidx<=4
            title(cornernames(spidx),'FontSize',11);
        end
        if ismember(spidx,[1:4:9])
            text(-500,1620/2,rcpnames(sspi),'HorizontalAlignment','center','rotation',90,'fontweight','bold','FontSize',11);
            sspi=sspi+1;
        end

        % add colorbar at the bottom
        if spidx==9
            cb= mycbar(sprintf("Change in %s potential in %s future (%%)",potnames_cl5{p},tf_full(seltf)),'southoutside');
            cb.Position=[0.09 0.06 0.9 0.015];    % shift it to bottom of fig
        end

        % add alphabetlabel
        if spidx==1
            addAlphaLabel(2,"outside")
        end
        c=c+1;
    end
end
