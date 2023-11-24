% Prepare cost plots with basin level for 3x2 energy (Tech+sust) and 2 hazard rep cases
% Prepare cost plots with sub-basin wise distinction of mixed sust case

clc
clearvars
close all
run('myVarNames.m')

load(fullfile(pwd,'\output\Figs_trial\MainScenarios.mat'), 'pcsout',...
    'runnames','pcsouthaz',  'runnameshaz','basindata','catchments_cl')
pcsout=[pcsout pcsouthaz];
runnames=[runnames runnameshaz];

%% FINAL: Loop FULL + Hazard Cost curve w points for DP and RP - assumes tech and sust scenarios are back to back
fig=figure('Position', get(0, 'Screensize'),'color','w');

s(1)=subplot(1,3,1:2);
runnames_cl=strrep(strrep(extractAfter(runnames,'Full_'),'_Tech_Fin','-Technical'),'_Sust_','-Sustainable-'); %,["Large-", "Medium-", "Mixed-"]);
myscens=[2 1 4 3 6 5 13 14];
copts=[];
%create my color palette with one shade lighter version for Tech case
for i=1:3 % 3 scens % create lighter shade for each scenario
    copts=[copts; cmap3_mainenergy(i,:);  brighten(cmap3_mainenergy(i,:),0.6);  ];
end

% Add greys for hazard
copts=[copts; flipud(cbrewer2( 'Greys', 4))];
c=1;
mylegend={};
mylegend_scname={};
for scn= myscens
    scol=copts(c,:);c=c+1; %bcoz i may not go serially
    plotCostCurve(pcsout{scn}.COEAlls, pcsout{scn}.PnetAlls, pcsout{scn}.SysIDs,scol,'-',.8,runnames_cl{scn});
    %    plotCostCurveScatter(pcsout{i}.COEAlls, pcsout{i}.PnetAlls, pcsout{i}.SysIDs,scol,1);
    mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
    mylegend_scname= [mylegend_scname(:)' {''} {''} runnames_cl(scn) ];
end
% Add lines for financial cutoffs
yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
text(250,0.03, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)

xlim([0,305])
ylabel('Per unit production cost (USD2010 per KWh)',FontWeight='bold')
xlabel('Basin-wide cumulative energy (TWh per year)',FontWeight='bold')
set(gca,'YScale','log')
set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
%set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
grid on; box on

%lgd=legend(mylegend,'NumColumns',2,'Orientation','horizontal','Location','best');
%lgd=legend(mylegend_scname,'NumColumns',3,'Orientation','vertical', 'Interpreter','none','Location','north'); %
lgdbox=[0.3,0.88,0.32239582768331,0.07];
lgd=legend([planttypes {'   '} repmat([{''} {''} {'   '}],1,3) {''} {''} "   Sustainable" {''} {''} "   Technical" {''} {''} "   Cost-Based" {''} {''} "   Multi-hazard" ],...
    'NumColumns',5,'Orientation','vertical','Location','northoutside', 'Interpreter','none',...
    'EdgeColor','none','Position',lgdbox);
annotation('rectangle', lgdbox, 'LineWidth',0.3);
title(lgd,sprintf("        PLANT TYPE          %s   %s   %s    SCENARIO            RISK TYPE",searchtypes{:}))
myylim=ylim;

%% FINAL: Sub-basinwise cost curves for FULL Mixed sust
%figure
s(2)=subplot(1,3,3);
cmap8=cmap8_wong;
scn=6; % the Full mixed sust scenario
addfill=0;
mylegend={};
mylegend_basname={};
c=1;
for basID= 101:108
    scol=cmap8(c,:);c=c+1; %bcoz i may not go serially
    selbas=pcsout{scn}.subbasinIDs==basID;
    plotCostCurveScatter(pcsout{scn}.COEAlls(selbas), pcsout{scn}.PnetAlls(selbas), pcsout{scn}.SysIDs(selbas),scol,addfill,'-',basindata.basinnames{basID-100});
    mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
    mylegend_basname= [mylegend_basname(:)' {''} {''} [basindata.basinnames(basID-100) ]];
end

% Add lines for financial cutoffs
yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
text(18,0.03, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)
xlabel('Subbasin-wise cumulative energy (TWh per year)',FontWeight='bold')
ylabel('Per unit production cost (USD2010 per KWh)',FontWeight='bold')
ylim(myylim)
set(gca,'YScale','log')
set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
grid on; box on
title("Mixed sust")

% subcatchment map overlaid
%subplot(4,1,1)
%    axes('Position', [1-.3-.095, .75, .2, .2]) % for independent fig
s(3)=axes('Position', [1-.277, .8, .15, .15]); % for subplot case
h=imagescnan(catchments_cl);
colormap(cmap8)
set(h, 'AlphaData', .8*(catchments_cl>0))

% add names for subbasins
text(bas_x(1:8),bas_y(1:8),basindata.basinnames,"HorizontalAlignment","center","FontAngle","italic")

%     cb=colorbar('Ticks', basindata.basinIDs,'TickLabels',basindata.basinnames,'Direction','reverse','Location','westoutside');
%     % Make the colorbar transparent - Need to do this manually
%     cdata = cb.Face.Texture.CData;
%     cdata(end,:) = uint8(myalpha * cdata(end,:));
%     cb.Face.Texture.ColorType = 'truecoloralpha';
%     cb.Face.Texture.CData = cdata;

%%
set(s(1),'Position', [0.13          0.11          0.495          0.85]);
set(s(2),'Position', [0.665 .11 .2 .85]);
set(s(3),'Position', [0.69          0.8          0.15          0.15]);

% For labels manually set y positions as follows
%Mixed=1
%Medium=10
%Large=60

%% FINAL: Sub-basinwise cost curves for FULL Mixed TECH
figure
cmap8=cmap8_wong;
scn=5; % the Full mixed sust scenario
addfill=0;
mylegend={};
mylegend_basname={};
c=1;
for basID= 101:108
    scol=cmap8(c,:);c=c+1; %bcoz i may not go serially
    selbas=pcsout{scn}.subbasinIDs==basID;
    plotCostCurveScatter(pcsout{scn}.COEAlls(selbas), pcsout{scn}.PnetAlls(selbas), pcsout{scn}.SysIDs(selbas),scol,addfill,'-',basindata.basinnames{basID-100});
    mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
    mylegend_basname= [mylegend_basname(:)' {''} {''} [basindata.basinnames(basID-100) ]];
end

% Add lines for financial cutoffs
yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
text(18,0.03, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)
xlabel('Subbasin-wise cumulative energy (TWh per year)',FontWeight='bold')
ylabel('Per unit production cost (USD2010 per KWh)',FontWeight='bold')
ylim(myylim)
set(gca,'YScale','log')
set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
grid on; box on
title("Mixed TECH")

% subcatchment map overlaid
%subplot(4,1,1)
%    axes('Position', [1-.3-.095, .75, .2, .2]) % for independent fig
s(3)=axes('Position', [1-.277, .8, .15, .15]); % for subplot case
h=imagescnan(catchments_cl);
colormap(cmap8)
set(h, 'AlphaData', .8*(catchments_cl>0))

% add names for subbasins
text(bas_x(1:8),bas_y(1:8),basindata.basinnames,"HorizontalAlignment","center","FontAngle","italic")

%     cb=colorbar('Ticks', basindata.basinIDs,'TickLabels',basindata.basinnames,'Direction','reverse','Location','westoutside');
%     % Make the colorbar transparent - Need to do this manually
%     cdata = cb.Face.Texture.CData;
%     cdata(end,:) = uint8(myalpha * cdata(end,:));
%     cb.Face.Texture.ColorType = 'truecoloralpha';
%     cb.Face.Texture.CData = cdata;


%% Create subbasin wise bubble plot for just plant sizes
cmap8=cmap8_wong;
scn=6; % the Full mixed sust scenario
c=1;
figure
subplot(1,3,1)
hold on
for basID= 101:108
    scol=cmap8(c,:);c=c+1; %bcoz id may not go serially
    selbas=pcsout{scn}.subbasinIDs==basID;
    yy=sort(pcsout{scn}.PnetAlls(selbas));
    %plot(sort(pcsout{scn}.PnetAlls(selbas)),'.','DisplayName',basindata.basinnames{basID-100});
    bubblechart(1:length(yy),yy,yy,scol,'DisplayName',basindata.basinnames{basID-100},'MarkerFaceAlpha',0.20);
end
xlabel("CellID")
ylabel("GWh per year")
grid on; box on
legend
blgd = bubblelegend('GWh/yr');

% Create subbasin wise projects box plot
hold on
yy4box=[];
xx4box=[];
subbas=[];
c=1;
for basID= 101:108
    scol=cmap8(c,:);c=c+1; %bcoz id may not go serially
    selbas=pcsout{scn}.subbasinIDs==basID;
    yy4box=[yy4box;sort(pcsout{scn}.PnetAlls(selbas))];
    xx4box=[xx4box;sort(pcsout{scn}.COEAlls(selbas))];

    subbas=[subbas; ones(length(sort(pcsout{scn}.PnetAlls(selbas))),1)*basID];
end

subplot(1,3,2)
boxplot(yy4box,subbas,'Symbol','.','Whisker',Inf,'label',basindata.basinnames)%,scol,'DisplayName',basindata.basinnames{basID-100},'MarkerFaceAlpha',0.20);
set(gca, 'YScale', 'log')
xlabel("Subbasin")
ylabel("GWh per year")
grid on; box on

subplot(1,3,3)
boxplot(xx4box,subbas,'Symbol','.','Whisker',Inf,'label',basindata.basinnames)%,scol,'DisplayName',basindata.basinnames{basID-100},'MarkerFaceAlpha',0.20);
xlabel("Subbasin")
ylabel("COE in USD2010/kWh")
grid on; box on
yline(0.1)

%% Create subbasin wise projects swarm chart or histogram- but this does not look that good so pass
%subplot(1,2,2)
figure
hold on
yy4box=[];
subbas=[];
c=1;
for basID= 101:108
    scol=cmap8(c,:);c=c+1; %bcoz id may not go serially
    selbas=pcsout{scn}.subbasinIDs==basID;
    %subplot(4,4,c);scatterhist(pcsout{scn}.PnetAlls(selbas),pcsout{scn}.COEAlls(selbas),'Kernel','on')
    yy=sort(pcsout{scn}.PnetAlls(selbas));
    subplot(2,4,c-1);histogram(yy,'DisplayName',basindata.basinnames{basID-100})
    %swarmchart(basID,yy,5,scol)
    legend
    xlabel("GWh per year")
    xlim([0 4000])
     ylim([0 1000])
    set(gca, 'YScale', 'log', 'XScale', 'log')
grid on
end

%xlabel("Subbasin")
grid on; box on
