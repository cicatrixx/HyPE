% Fig_BasinTotalsBar
% Created By    : Sanita Dhaubanjar
% Created For	: SustainIndus project WP2 Paper publication
%=========================
% Code for loading pre-processed output files from Hydrus model to create
% final plots used for publication

clc
clearvars
close all
addpath(genpath(fullfile(pwd,"Hydrus")))
% Add folder path
run('myVarNames_Fut.m')
rootofigdata=fullfile(rootof,"FutFigs_cl0","FutFigs_cl1","FigData");

charlbl =  compose("(%s)",('a':'f').');

%% Load pre-proceesed data
% Load hist and future tech, fin, sust potential
load(fullfile(rootofigdata, 'Compile_FutScenario_Totals00.mat'),'tot_allscen','histHP','potlist')
futHP_allscen=tot_allscen;
histHP_mixed=histHP;

% Load hist and future theoretical potential
load(fullfile(rootofigdata,'FutTheoreticalPot_subbasintotals.mat'), 'subPot_GWh')
theorypot_subbasin_GWh=subPot_GWh; % last row in this is total of all basins, other rows are subbasin vals
disp("Loaded data")
clear tot_allscen histHP subPot_GWh

%% FINAL: Vertical bar plot w total and RP/DP and hist bar
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.005],[.08 .08],[.08 .08]); %[vgap hgap], marg_h -[low up],marg_w -[l r]
%Matlab doesnt allow double axis plots to be rotated w the cam option so manually rotate elsewhere
% Add zero cell to create gaps
theorypot_subbasin_GWh(9,26)=0;
sel2plot=[1,26,2:13,26,14:25];

figure
colororder(c_ssp26(sel2plot,:))
% First plot theoretical
plabel=1;
subplot(6,1,6)
% plot mid and far terms as separate
bar(diag(theorypot_subbasin_GWh(9,sel2plot)/1000'),'stacked','EdgeAlpha',0)
yline(theorypot_subbasin_GWh(9,1)/1000,'color','black','LineWidth',1.5,'LineStyle',':')
xline([2,15],'color',mygraylines,'LineWidth',1.05) % tf splits
set(gca,'YDir','reverse','XTick',1:25,'XTickLabels',[]) %,'YColor',[0.00,0.45,0.74]...
%'XAxisLocation','top') %,'XLim',[0 110]),'YDir','reverse',
ylabel(strcat(charlbl(plabel)," ",pottypes5(1)),'FontWeight','bold') % "TWh per yr"])
ylim([0 2500])
applymyplotformat('')

% Then plot other pots
pp=4;
for pidx=[1,2,4]
    plabel=plabel+1;
    pt=newpots4(pidx);
    selp=find(futHP_allscen.pottype==pt);

    % add gap
    selrows=[97; selp(1:12); 97; selp(13:24)];

    % plot total for hist and future
    subplot(6,1,pp)
    idata=[histHP_mixed.("Total TWh/yr")(histHP_mixed.pottype==pt);   futHP_allscen.("Total TWh/yr")(selrows)]; % add historical
    b1=bar(diag(idata),'stacked','EdgeAlpha',0);
    hold all
    % add boxes for RP/DP totals for hist and future
    idata=[histHP_mixed{(histHP_mixed.pottype==pt), 2:3};   futHP_allscen{selrows, 2:3}]; % add historical
    b2=bar(idata,.7,  'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.5,'LineStyle','-');
    ylabel(strcat(charlbl(plabel)," ",pottypes3(min(pidx,3))),'FontWeight','bold' )

    % Change color of  bar
    for k = 1:2
        b2(k).EdgeColor = mycolors{k}*0.6;
    end

    % add lines for hist pot and separaters for tfs
    b3=yline(futHP_allscen.("Historical TWh/yr")(selp(1)),'color','black','LineWidth',1.5,'LineStyle',':');
    xticks(1:(length(selrows)+1))
    xline([2,15],'color',mygraylines,'LineWidth',1.05) % tf splits
    ylim([0 440])

    % add diamonds for number of RP/DP hist and then fut
    yyaxis right
    scshift=0;
    idata=[histHP_mixed{(histHP_mixed.pottype==pt), 5}; futHP_allscen{selrows, 5}];
    s1(1)=scatter([1:length(selrows)+1]-scshift,idata,30,cl_RP*.6,'d','filled');
    idata=[histHP_mixed{(histHP_mixed.pottype==pt), 5}+histHP_mixed{(histHP_mixed.pottype==pt), 6}; futHP_allscen{selrows, 6}+futHP_allscen{selrows, 5}];
    s1(2)=scatter([1:length(selrows)+1]-scshift,idata,30, cl_DP*.6,'d','filled'); %manually stack scatter

    % add unit label on the left and right
    if pidx==2
        yyaxis left
        text(-1,210,"Energy in TWh per year",'HorizontalAlignment','center','FontWeight','bold','Rotation',90)
        yyaxis right
        ylabel('Number of projects','FontWeight','bold')

    end
    ax = gca;
    ax.YAxis(2).Color = extraaxis ;

    % add uniform  yyaxis only for fin and sust cos it is hard to see diff
    % between RP and DP in them if i use same axis as for tech
    if pp<=3
        ylim([0 3000])
    end
    applymyplotformat('')

    % rotate subplot
    % camroll(-90)

    % add xticklabels
    if pp==4
        xticklabels([histrcpcornernames(1) '' cornernamesrep '' cornernamesrep])
        xtickangle(90)
        % text([7.5 21.5],-400*[1 1],repmat(strjoin(rcpnames,'         '),1,2),'HorizontalAlignment','center','FontSize',11)
        yyaxis left
        % add rcp labels
        text([2+[2.5:4:12] 15+[2.5:4:12]],-300*[1 1 1 1 1 1],[rcpnames rcpnames],'HorizontalAlignment','center','FontSize',10,'FontWeight','bold')
        % add mid far labels
        text([8.5 15+6.5],-400*[1 1],strcat(tname,": ", tframe),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold')%,'Rotation',90)

    else
        xticklabels([])
    end

    pp=pp-1;
end
%add legend
l=legend([b1(3:14) b2 s1 b1(1)],[repelem({'.',},12-4) cornernames...
    'RP energy' 'DP energy' 'RP number' 'DP number' 'Historical'],'NumColumns',5, 'Location','northoutside');
l.Title.String=strjoin(rcpnames,'      ');
% rotate subplot so plot becomes horizontal bar
%camroll(-90)


%%
disp("***************************************EOF***************************************")

