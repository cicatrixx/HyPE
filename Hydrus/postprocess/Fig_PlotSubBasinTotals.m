%% Combine sub-basin wise onion and bar plots into one
clear all
close all
addpath(genpath('G:\SurfDrive\GitConnect\Hydrus\devFiles\'), ...
    genpath('G:\SurfDrive\HPmodel\Hydrus\'))
run('myVarNames.m')
load('G:\SurfDrive\HPmodel\output\Figs_trial\MainScenarios.mat', 'basindata','catchments_cl',...
    'subPot_GWh', 'subPot_num', 'subPot_GWh_perarea','subPot_GWh_percapita')
load('G:\SurfDrive\HPmodel\data\ASIA\Basin_UIB\PantpeBasin_103.mat', 'outside')
catchments_cl=maskBasin(catchments_cl,~outside);

cmap8=brighten(cmap8_wong,0);
label_subbasin_basin =[ basindata.basinnames{:},"ALL Indus"];

% %% Manually assign xy for onions - already in myVarnames
% %'Kabul'  'Swat'  'Indus''Jhelum''Chenab'    {'Ravi'}    {'Beas'}  {'Satluj'} 'All basins'
% bas_x =[670     1050     1800    1360   1600    1400    1900    2300 2400];
% bas_y =[410     380     500     600     890    1200    1100    1250 300] ;

%% bar Plot all Tech/econ/sust pot types as bars, add basinwide and existing as empty bars
selbasins=1:9;
figure
%    cmap5_orange=cbrewer2('Oranges',5); %lines;
selpots=2:4; %size(subPot_GWh,2);
for sc_idx=1:3
    subplot(2,3,sc_idx)
    idata=subPot_GWh_perarea*1e6;
    %Add TFS
    b1=bar(idata(selbasins,selpots,sc_idx),'FaceAlpha',baralpha);
    if sc_idx==1;   ylabel({'Specific potential'; '(GWh/yr per km^2)'},'FontWeight','bold'); end
    xtickangle(30)
    xticklabels(label_subbasin_basin)
    applymyplotformat('',cmap6_pottype(selpots,:))
    set(gca,'XGrid','off')
    ylim([0, 2.6])
    hold on
    %Add visualized
    b2=bar(idata(selbasins,5,sc_idx),'LineStyle','-','FaceColor',cmap6_pottype(5,:),'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'DisplayName','Visualized');
    %Add existing
    b3=bar(idata(selbasins,6,sc_idx),'LineStyle','-','FaceColor',cmap6_pottype(6,:),'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'DisplayName','Existing');
    %b4=bar(repmat(idata(9,selpots,sc_idx),8,1),'LineStyle','-','FaceColor','k','EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'DisplayName','Basin level value');


    subplot(2,3,sc_idx+3)
    idata=subPot_GWh_percapita*1e3; % subPot_GWh_percapita*1e3; %
    %Add TFS
    b11=bar(idata(selbasins,selpots,sc_idx),'FaceAlpha',baralpha);
    if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
    xtickangle(30)
    xticklabels(label_subbasin_basin)
    applymyplotformat('',cmap6_pottype(selpots,:))
    set(gca,'XGrid','off')
    ylim([0, 33])
    hold on
    %Add visualized
    b21=bar(idata(selbasins,5,sc_idx),'LineStyle','-','FaceColor',cmap6_pottype(5,:),'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'DisplayName','Visualized');
    %Add existing
    b31=bar(idata(selbasins,6,sc_idx),'LineStyle','-','FaceColor',cmap6_pottype(6,:),'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'DisplayName','Existing');
    %b41=bar(repmat(idata(9,selpots,sc_idx),8,1),'LineStyle','-','FaceColor','k','EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'DisplayName','Basin level value');

    % Change color of  bar edges
    for k = 1:length(selpots)
        b1(k).EdgeColor = cmap6_pottype(selpots(k),:);
        b11(k).EdgeColor = cmap6_pottype(selpots(k),:);
    end
end
%    legend(pottypes6(selpots),'Orientation','horizontal','Location','north') %,'NumColumns',4  )
legend([pottypes6{[selpots 5 6]}, "All basin value"],'Orientation','horizontal','Location','north') %,'NumColumns',4  )

%% FINAL: Combine onion and bar
figure
% subplot (4, 3,1:6)
subtightplot (4, 3,1:6, [0.025 0.0],[.08 .02],[.08 .01]); %[vgap hgap], marg_h -[low up],marg_w -[l r]
selsearch=3;
my_bubble_size=subPot_GWh_percapita*1e3;
selcols=1:4;
h=imagescnan(catchments_cl);hold all
colormap(cmap8)
set(h, 'AlphaData', subbasinalpha*(catchments_cl>0))
%Add theoretical
bubblechart(bas_x,bas_y,my_bubble_size(:,1,selsearch),cmap6_pottype(1,:),"MarkerFaceAlpha",0.6)
%Add TFS
for ii=selcols(2:end)
    hold on
    bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0.3)
end
ii=ii+1;
%Add visualized
bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0.3)% ,'LineWidth',1.5)
%Add existing
ii=ii+1;
bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0.3)%,'LineWidth',1.5,'DisplayName',"Existing")
legend(pottypes6,'Location','northwest')
xticklabels(basindata.basinnames)
bubblesize([5, 120])
blgd=bubblelegend("MWh/yr per capita",'Style','telescopic','Box','off','Location','southwest');
title(sprintf("Sub-basin MWh/yr per capita for %s search",searchtypes{selsearch}))
text(bas_x,bas_y,label_subbasin_basin,'HorizontalAlignment','center','FontAngle','italic')

selbasins=[9 1:8];
%    cmap5_orange=cbrewer2('Oranges',5); %lines;
selpots=2:4; %size(subPot_GWh,2);
for sc_idx=1:3
    subtightplot(4,3,6+sc_idx, [0.02 0.01],[.08 .02],[.08 .01]); %[vgap hgap], marg_h -[low up],marg_w -[l r]
    idata=subPot_GWh_perarea*1e6;
    %Add TFS
    b1=bar(idata(selbasins,selpots,sc_idx),'FaceAlpha',baralpha);
    yticks([0:0.5:2.5])  ;   ylim([0, 2.8])
    if sc_idx==1;   ylabel({'Specific potential'; '(GWh/yr per km^2)'},'FontWeight','bold'); end
    if sc_idx>1;     yticklabels("");     end
    %xtickangle(30)
    xticklabels([])
    %xticklabels(label_subbasin_basin)
    applymyplotformat(searchtypes{sc_idx},cmap6_pottype(selpots,:))
    set(gca,'XGrid','off')
   
    hold on
    %Add visualized
    b2=bar(idata(selbasins,5,sc_idx),'LineStyle','-','FaceColor',cmap6_pottype(5,:),'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'DisplayName','Visualized');
    %Add existing
    b3=bar(idata(selbasins,6,sc_idx),'LineStyle','-','FaceColor',cmap6_pottype(6,:),'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'DisplayName','Existing');
    %b4=bar(repmat(idata(9,selpots,sc_idx),8,1),'LineStyle','-','FaceColor','k','EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'DisplayName','Basin level value');
    xline(1.5,'Color','k')

    subtightplot(4,3,6+sc_idx+3, [0.02 0.01],[.08 .02],[.08 .01]); %[vgap hgap], marg_h -[low up],marg_w -[l r]
    idata=subPot_GWh_percapita*1e3; % subPot_GWh_percapita*1e3; %
    %Add TFS
    b11=bar(idata(selbasins,selpots,sc_idx),'FaceAlpha',baralpha);
    yticks([0:5:30])  ;
    ylim([0, 33])
    if sc_idx==1;   ylabel({'Per capita potential'; '(MWh/yr per capita)'},'FontWeight','bold');end
    if sc_idx>1;       yticklabels("");     end
    %xtickangle(30)
    xticklabels(label_subbasin_basin(selbasins))
    applymyplotformat('',cmap6_pottype(selpots,:))
    set(gca,'XGrid','off')
    
    hold on
    %Add visualized
    b21=bar(idata(selbasins,5,sc_idx),'LineStyle','-','FaceColor',cmap6_pottype(5,:),'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'DisplayName','Visualized');
    %Add existing
    b31=bar(idata(selbasins,6,sc_idx),'LineStyle','-','FaceColor',cmap6_pottype(6,:),'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'DisplayName','Existing');
    %b41=bar(repmat(idata(9,selpots,sc_idx),8,1),'LineStyle','-','FaceColor','k','EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1,'DisplayName','Basin level value');
    xline(1.5,'Color','k')

    % Change color of  bar edges
    for k = 1:length(selpots)
        b1(k).EdgeColor = cmap6_pottype(selpots(k),:);
        b11(k).EdgeColor = cmap6_pottype(selpots(k),:);
    end
end
%legend(pottypes6{[selpots 5 6]},'Orientation','horizontal','Location','north') %,'NumColumns',4  )
