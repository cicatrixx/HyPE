% Fig_SubBasinPotentials_PerCapita
% Created By    : Sanita Dhaubanjar
% Created For	: SustainIndus project WP2 Paper publication
%=========================
% Code for loading pre-processed output files from Hydrus model to create
% final plots used for publication
% Plot sub-basin wise total potential and % change in potential for far
% future

clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

%% Load saved subbasin data for fut potentials
ofmatname=fullfile(rootof,"FutFigs_cl0","FutFigs_cl1","FigData",'MainFutScenarios.mat');
load(ofmatname,'potnames','catchments_cl', 'subPot_MWh_percapita_futpop','subPopMWhPerCapitaPrctChange')
disp("Loaded saved subPot data")
potnames_cl5=[potnames(1) newpots4cl];

%skip the sust-tech pot in excel
selpot=[1 2 3 5 6 7 ];

%% FINAL: Bar plot of total + spatial of far future for subbasin-wise potential - FOR ONE POTENTIAL ALL FUTS
%selscendata
for p=selpot(4) % 4=sustainable potential
    idata=squeeze(subPot_MWh_percapita_futpop(:,p,:)); %2 and 4
    idataprct=squeeze(subPopMWhPerCapitaPrctChange(:,p,:)); % basins are rowwise and scens are column wise - mid fut first, far fut second

    %barplot is easier to read than scatter. so kept bar
    figure;
    subplottight(5,4,1:4);
    b1=bar(idata(1:8,2:end),'EdgeAlpha',0);
    hold on
    b2=bar(idata(1:8,1),'FaceAlpha',0);
    xline(1:8,'color',mygraylines,'LineWidth',1.05) % tf splits
    xticklabels(basinlabels.basinnames)
    xtickangle(0)
    xlim([0.4,8.6])
    ylabel("MWh per capita per year",'fontweight','bold','FontSize',11)
    applymyplotformat(sprintf("Total %s potential",potnames_cl5{p}),c_ssp)

    % Add min energy threshold
    b3=yline(energyreq_MWhpercapita,':r');

    l=legend([b1(13:24), b2, b3],[repelem({' ',},24/2-4) cornernames 'Historical' 'Minimum energy threshold'],'NumColumns',4);
    l.Title.String=strjoin(rcpnames,'      ');

    % Add annotations
    addAlphaLabel(1,"outside")

    % add mid-far label
    topy=ylim;%text(1,max(idataTWh(1,:)),strcat("\leftarrow ", tname(1),"   " ,tname(2)," \rightarrow "),'HorizontalAlignment','center')%,'FontAngle','italic')
    text(1,floor(topy(2)),strcat(tname(1),"      " ,tname(2)),'HorizontalAlignment','center')%,'FontAngle','italic')

    % Shift downwards and increase height a bit
    h=gca;
    h.Position=[h.Position(1) h.Position(2)-.05 h.Position(3) h.Position(4)+.05];

    % Add spatial plots for far future
    seltf=2;
    %get min-max for far future scenarios to set limits for all cmaps
    mylim=[floor(min(idataprct(1:8,14:25),[],'all')) ceil(max(idataprct(1:8,14:25),[],'all'))]

    if p==selpot(4) % For sustainable pot
        mylim=[-79 106];
        %       colormap(cbrewer2('RdYlBu'))
        tmpcmap=[brighten(flipud(cbrewer2('YlOrRd',abs(mylim(1))+10)),0.1);  cbrewer2('Blues',abs(mylim(2))+5)];
        colormap(tmpcmap([10:abs(mylim(1))+10,abs(mylim(1))+11:end],:))
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
            text(-500,1620/2,sspnames(sspi),'HorizontalAlignment','center','rotation',90,'fontweight','bold','FontSize',11);
            sspi=sspi+1;
        end

        % add colorbar at the bottom
        if spidx==9
            cb=mycbar(sprintf("%% change in MWh per capita per year for %s potential in %s future",potnames_cl5{p},tf_full(seltf)),'southoutside');
            cb.Position=[0.09 0.06 0.9 0.015];    % shift it to bottom of fig
        end

        % add alphabetlabel
        if spidx==1
            addAlphaLabel(2,"outside")
        end
        c=c+1;
    end
end
