% Prepare data for cost plots - try out different types

% Load data for main scenarios - same as gettotals
clc
clearvars
close all
addpath(genpath('G:\SurfDrive\GitConnect\Hydrus\devFiles\'), ...
    genpath('G:\SurfDrive\HPmodel\Hydrus\'))

run('myVarNames.m')

load('G:\SurfDrive\HPmodel\output\Figs_trial\MainScenarios.mat', 'pcsout',...
    'runnames','pcsouthaz',  'runnameshaz','basindata','catchments_cl')
 pcsout=[pcsout pcsouthaz];
runnames=[runnames runnameshaz];

ifname = sprintf('%s\\data\\%s\\Basin_UIB\\PantpeBasin_101.mat', rootf, continent_in);
datain=load(ifname,'outside','channel_main_trib');

%% GOOD: Old layout Loop to make 2 figs, mult subplot of dist of size + Color coded size vs cost vs type
myscens=[1:6];
nplot=length(myscens);

%subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.03]);
fig_order= [1 1 2 2 3 3 4 4 5 5];      %    [1,3,5,2,4,6]; %

for i=1:length(myscens)
    mydata=pcsout{i};
    %mydata.res=pcsin.res;
    % Plot figs w 3 panels
    if rem(i,2)
        figure(10) %for tech
    else
        figure(20) %for sust
    end
    
    subplot(1,nplot/2,fig_order(i))
    plot(mydata.COEAlls(mydata.SysIDs==1),mydata.PnetAlls(mydata.SysIDs==1),'.','Color',cl_DP,'DisplayName',planttypes{2})
    hold all
    plot( mydata.COEAlls(mydata.SysIDs==2),mydata.PnetAlls(mydata.SysIDs==2),'.','Color',cl_RP,'DisplayName',planttypes{1})
    %title("Cost of electricity for various projects")
    
    %add lines for different categories of HP
    yline(HPsz_greaterthan(1),'--k');
    yline(HPsz_greaterthan(2),'--k');
    yline(HPsz_greaterthan(3),'--k');
    yline(HPsz_greaterthan(4),'--k');
    yline(HPsz_greaterthan(5),'--k');
    %xline(Techlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Technical potential limit'); %Tech potential
    %text(Techlim*1.15, 0,sprintf("\\it Technical potential <= %0.2f USD/kWh",Techlim),'FontSize',12)
    
    xline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(costlim*1.15, 0, sprintf("\\it Economic potential\n <= %0.2f USD/kWh",costlim),'FontSize',12)
        set(gca, 'YScale', 'log','XScale', 'log')

    ylabel('GWh','FontSize',10)
    xlabel('Per unit production cost (USD2010 per KWh)','FontSize',10)
    grid on
    %title(char(64+i),'FontName','franklin gothic medium','FontSize',14)
    title(runnames{myscens(i)},'Interpreter', 'none')
    ylim([0.02 4*10^4])
    %title(sprintf("For %s in Upper Indus",selrun),'Interpreter', 'none' )
    %sgtitle(sprintf("For %s in Upper Indus",selrun),'Interpreter', 'none' )
end
%add labels for different categories of HP
figure(10);text(repmat(max(mydata.COEAlls),length(HPclass),1),[HPsz_greaterthan*5 HPsz_greaterthan(end)*.9] ,HPclass,'FontWeight','bold')
legend('Diversion type','River Power')

figure(20);text(repmat(max(mydata.COEAlls),length(HPclass),1),[HPsz_greaterthan*5 HPsz_greaterthan(end)*.9] ,HPclass,'FontWeight','bold')
legend('Diversion type','River Power')



%% Try scatterhist standard - cant get log axis to work in the subplots :E
i=5;
mydata=pcsout{i};
mydata.SysIDs=categorical(mydata.SysIDs,[2 1], planttypes);
h=scatterhist(mydata.COEAlls,mydata.PnetAlls,'Group',mydata.SysIDs,'Kernel','on','Location','SouthEast',...
    'Direction','out','Color','bb','LineStyle',{':','--'},...
    'LineWidth',[1,1,1],'Marker','+od','MarkerSize',[4,5,6]);
 grid minor
%set(h(1), 'YScale', 'log','XScale', 'log')
% 2= bottom, 3= right
%set(h(2:3),'XScale', 'log');

%% GOOD: Try scatterhist w box plots from mlab example-w scen and ptype grouping - One scenario!
i=5;
mydata=pcsout{i};
mydata.SysIDs=categorical(mydata.SysIDs,[2 1], planttypes);
grpscenptyp = strrep(strcat(extractAfter(mydata.runname,'Energy_'),"-",string(mydata.SysIDs)),'_','-');

figure
h=scatterhist(mydata.COEAlls,mydata.PnetAlls,'Group',grpscenptyp,'Kernel','on','Location','SouthEast',...
    'Direction','out','Color',cmap3energy,... %'LineStyle',{':','--'},'LineWidth',[1,1,1],...
    'Marker','o^','MarkerSize',6);
 grid minor
set(h(1), 'YScale', 'log','XScale', 'log')
% 2= bottom, 3= right
%set(h(2:3),'XScale', 'log');
%box(h(2:3),'on')
% add boxplot to scatter
boxplot(h(2),mydata.COEAlls,mydata.SysIDs,'orientation','horizontal',...
     'label',{'',''},'color',cmap3energy,'PlotStyle','compact','OutlierSize',2,'Symbol','.');

boxplot(h(3),mydata.PnetAlls,mydata.SysIDs,'orientation','horizontal',...
     'label', {'',''},'color',cmap3energy,'PlotStyle','compact','OutlierSize',2,'Symbol','.');
view(h(3),[270,90]);  % Rotate the Y plot on the E

set(h(2:3),'XScale','log','XTickLabel','')
grid(h(2:3),'minor')

% Add annotations
%addHPcategoryLines()
xline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
text(costlim*1.15, 0, sprintf("\\it Economic potential\n <= %0.2f USD/kWh",costlim),'FontSize',12)

ylabel('GWh','FontSize',10)
xlabel('Per unit production cost (USD2010 per KWh)','FontSize',10)
grid on


%% GOOD: Try scatterhist w box plots-w scen and ptype grouping - add 3 energy tech econs!
myscens=[1, 3,5];
COEAllsComp=[];
PnetAllsComp=[];
SysIDsComp=[];
GrpscenptypComp=[];
for i=myscens
    mydata=pcsout{i};
    mydata.SysIDs=categorical(mydata.SysIDs,[2 1], planttypes);
    grpscenptyp = erase(strrep(strcat(extractAfter(mydata.runname,'Energy_'),"   |  ",string(mydata.SysIDs)),'_','-'),["Full-", "-Tech-Fin"]);
    COEAllsComp =[COEAllsComp;mydata.COEAlls];
    PnetAllsComp=[PnetAllsComp;mydata.PnetAlls];
    SysIDsComp=[SysIDsComp;mydata.SysIDs];
    GrpscenptypComp=[GrpscenptypComp; grpscenptyp];
end


figure
h=scatterhist(COEAllsComp,PnetAllsComp,'Group',GrpscenptypComp,'Kernel','on','Location','SouthEast',...
    'Direction','out','Color',cmap3energy,... %'LineStyle',{':','--'},'LineWidth',[1,1,1],...
    'Marker','o^','MarkerSize',6);
 grid minor
set(h(1), 'YScale', 'log','XScale', 'log')
% 2= bottom, 3= right
%set(h(2:3),'XScale', 'log');
%box(h(2:3),'on')
% add boxplot to scatter
boxplot(h(2),COEAllsComp,GrpscenptypComp,'orientation','horizontal',...
     'color',cmap3energy,'PlotStyle','compact','OutlierSize',2,'Symbol','.'); %'label',{'',''},

boxplot(h(3),PnetAllsComp,GrpscenptypComp,'orientation','horizontal',...
     'color',cmap3energy,'PlotStyle','compact','OutlierSize',2,'Symbol','.');
view(h(3),[270,90]);  % Rotate the Y plot on the E

set(h(2:3),'XScale','log','XTickLabel','')
grid(h(2:3),'minor')

% Add annotations
%addHPcategoryLines()
xline(costlim,'LineStyle','--','Color','b','LineWidth',1.25,'DisplayName','Financial potential limit'); %Economic potential
%text(costlim*1.15, 0, sprintf("\\it Economic potential\n <= %0.2f USD/kWh",costlim),'FontSize',12)
% xline(h(2),costlim,'LineStyle','--','Color','b','LineWidth',1.25,'DisplayName',''); %Economic potential
ylabel('GWh','FontWeight','bold')
xlabel('Per unit production cost (USD2010 per KWh)','FontWeight','bold')
grid on

l=legend(); 
l.Title.String="SCENARIO             PLANT TYPE";
sgtitle("Full Tech Econ Scenarios")

%% GOOD: Try scatterhist w box plots-w scen and ptype grouping - add 3 energy susts!
myscens=[2,4,6];
COEAllsComp=[];
PnetAllsComp=[];
SysIDsComp=[];
GrpscenptypComp=[];
for i=myscens
    mydata=pcsout{i};
    mydata.SysIDs=categorical(mydata.SysIDs,[2 1], planttypes);
    grpscenptyp = erase(strrep(strcat(extractAfter(mydata.runname,'Energy_'),"   |  ",string(mydata.SysIDs)),'_','-'),["Full-", "-Tech-Fin","-Sust-RiskAverse"]);
    COEAllsComp =[COEAllsComp;mydata.COEAlls];
    PnetAllsComp=[PnetAllsComp;mydata.PnetAlls];
    SysIDsComp=[SysIDsComp;mydata.SysIDs];
    GrpscenptypComp=[GrpscenptypComp; grpscenptyp];
end


figure
h=scatterhist(COEAllsComp,PnetAllsComp,'Group',GrpscenptypComp,'Kernel','on','Location','SouthEast',...
    'Direction','out','Color',cmap3energy,... %'LineStyle',{':','--'},'LineWidth',[1,1,1],...
    'Marker','o^','MarkerSize',6);
 grid minor
set(h(1), 'YScale', 'log','XScale', 'log')
% 2= bottom, 3= right
%set(h(2:3),'XScale', 'log');
%box(h(2:3),'on')
% add boxplot to scatter
boxplot(h(2),COEAllsComp,GrpscenptypComp,'orientation','horizontal',...
     'color',cmap3energy,'PlotStyle','compact','OutlierSize',2,'Symbol','.'); %'label',{'',''},

boxplot(h(3),PnetAllsComp,GrpscenptypComp,'orientation','horizontal',...
     'color',cmap3energy,'PlotStyle','compact','OutlierSize',2,'Symbol','.');
view(h(3),[270,90]);  % Rotate the Y plot on the E

set(h(2:3),'XScale','log') %,'XTickLabel','')
grid(h(2:3),'minor')

% Add annotations
%addHPcategoryLines()
xline(costlim,'LineStyle','--','Color','b','LineWidth',1.25,'DisplayName','Financial potential limit'); %Economic potential
%text(costlim*1.15, 0, sprintf("\\it Economic potential\n <= %0.2f USD/kWh",costlim),'FontSize',12)
% xline(h(2),costlim,'LineStyle','--','Color','b','LineWidth',1.25,'DisplayName',''); %Economic potential
ylabel('GWh','FontWeight','bold')
xlabel('Per unit production cost (USD2010 per KWh)','FontWeight','bold')
grid on

l=legend(); 
l.Title.String="SCENARIO             PLANT TYPE";
sgtitle("Full Sust Scenarios")


%% Try scatterhistogram for one scenario- bad bcos not customizable even no log axis
validvars=~isnan(mydata.COEAlls);
tbl = table(repmat(mydata.runname,sum(validvars),1),mydata.COEAlls(validvars),mydata.PnetAlls(validvars),mydata.SysIDs(validvars),...
    'VariableNames',{ 'Scenario','COE [USD2010/kWh]', 'PnetAlls [GWh/yr]','Plant type'});

[idx,scentype,sptype] = findgroups(string(tbl.Scenario),string(tbl.("Plant type")));
tbl.grpscenptyp = strcat(scentype(idx),"-",sptype(idx));
figure
s = scatterhistogram(tbl,'COE [USD2010/kWh]', 'PnetAlls [GWh/yr]', ...
    'GroupVariable','grpscenptyp','HistogramDisplayStyle','smooth', ...
    'LineStyle','-');

%% Try gscatter - doesnt work here but otherwise is a promising use of color and symbols if not needing fancy format controls
x1=[];
x=cbrewer2('set3',8);
for i=1:8
    x1=[x1; repmat(x(i,:),2,1)];
end
figure
%pcsout{i}.COEAlls(selbas), pcsout{i}.PnetAlls(selbas), pcsout{i}.SysIDs(selbas)
gscatter(pcsout{i}.COEAlls, pcsout{i}.PnetAlls, pcsout{i}.SysIDs,...
     grpscenptyp, x1 ,25)
 

%% Create histogram plot
 grps=unique(GrpscenptypComp(~ismissing(GrpscenptypComp)));
 figure
 hold all
 for i=1:length(grps)
    histogram(COEAllsComp(GrpscenptypComp==grps(i)))
 end
%% Create violin like swarm plot using xy data
grps=unique(GrpscenptypComp(~ismissing(GrpscenptypComp)));
figure
swarmchart(categorical(GrpscenptypComp,grps),PnetAllsComp,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
set(gca,'YScale', 'log');
box on
grid on

%% GOOD: Create violin like swarm plot using table  - 3 energy Sust 
myscens=[2,4,6];
COEAllsComp=[];
PnetAllsComp=[];
SysIDsComp=[];
GrpscenptypComp=[];
for i=myscens
    mydata=pcsout{i};
    mydata.SysIDs=categorical(mydata.SysIDs,[2 1], planttypes);
    grpscenptyp = erase(strrep(strcat(extractAfter(mydata.runname,'Energy_'),"   |  ",string(mydata.SysIDs)),'_','-'),["Full-", "-Tech-Fin","-Sust-RiskAverse"]);
    COEAllsComp =[COEAllsComp;mydata.COEAlls];
    PnetAllsComp=[PnetAllsComp;mydata.PnetAlls];
    SysIDsComp=[SysIDsComp;mydata.SysIDs];
    GrpscenptypComp=[GrpscenptypComp; grpscenptyp];
end

tbl_Sust = table(categorical(GrpscenptypComp,grps),COEAllsComp,PnetAllsComp,...
    'VariableNames',{ 'Scenario-plant','COE [USD2010/kWh]', 'PnetAlls [GWh/yr]'});

figure
%tiledlayout('flow')
%nexttile;
s1 = swarmchart(tbl_Sust,'Scenario-plant','PnetAlls [GWh/yr]','ColorVariable','COE [USD2010/kWh]','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
set(gca,'YScale', 'log');
mycbar("COE [USD2010/kWh]")
box on
grid on
s1.SizeData = 5;
addHPcategoryLines(1,'h')

% same color coded by scenarios
%nexttile; cant do same fig as color map can be different!
figure
s2 = swarmchart(tbl_Sust,'Scenario-plant','PnetAlls [GWh/yr]','ColorVariable','Scenario-plant','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
set(gca,'YScale', 'log','ColorScale', 'log');
s2.SizeData = 5;
box on
grid on
colormap(cmap3energy)
addHPcategoryLines(1,'h')

%% Create violin swarm + Box - cant seem to overlay them because the axis is considered different
figure
subplot(2,1,1)
s2 = swarmchart(tbl_Sust,'Scenario-plant','PnetAlls [GWh/yr]','ColorVariable','Scenario-plant','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
set(gca,'YScale', 'log');
s2.SizeData = 5;
box on
grid on
colormap(cmap3energy)
%addHPcategoryLines(1,'h')
subplot(2,1,2)
boxplot(tbl_Sust{:,"PnetAlls [GWh/yr]"},tbl_Sust{:,"Scenario-plant"},'label',{'','','','','','',},...%'orientation','horizontal',...
     'color',cmap3energy,'PlotStyle','compact','OutlierSize',2,'Symbol','','FactorDirection','list','whisker', inf); 
set(gca,'YScale', 'log');
hold on
grid on

%swarmchart(tbl_Sust,"Scenario-plant","PnetAlls [GWh/yr]")

%% Create color for 3 energy scenarios
cmap3energy=[];
x=cbrewer2('Set2',3);  % dark2 or set2 are color blind and print friendly
for i=1:3 % 3 scens
    cmap3energy=[cmap3energy; x(i,:); brighten(x(i,:),-0.6) ];
end

%% Create table of 3 energy Sust 
myscens=[2,4,6];
ScensComp=[];
COEAllsComp=[];
PnetAllsComp=[];
SysIDsComp=[];
GrpscenptypComp=[];
for i=myscens
    mydata=pcsout{i};
    mydata.SysIDs=categorical(mydata.SysIDs,[2 1], planttypes);
    ndata=numel(mydata.SysIDs);
    scenname=erase(strrep(extractAfter(mydata.runname,'Energy_'),'_','-'),["Full-", "-Tech-Fin","-Sust-RiskAverse"]);
    grpscenptyp = strcat(scenname ,"|",string(mydata.SysIDs));
    ScensComp=[ScensComp; extractBefore(grpscenptyp,"|")];
    COEAllsComp =[COEAllsComp;mydata.COEAlls];
    PnetAllsComp=[PnetAllsComp;mydata.PnetAlls];
    SysIDsComp=[SysIDsComp;mydata.SysIDs];
    GrpscenptypComp=[GrpscenptypComp; grpscenptyp];
end

 grps=unique(GrpscenptypComp(~ismissing(GrpscenptypComp)));
grps=grps([2 1 4 3 6 5]);
tbl_Sust_all = table(categorical(GrpscenptypComp,grps),COEAllsComp,PnetAllsComp,SysIDsComp,categorical(ScensComp,searchtypes),...
    'VariableNames',{ 'Scenario-plant','COE [USD2010/kWh]', 'PnetAlls [GWh/yr]','Ptype','Scen'});
%select only not nans
tbl_Sust=tbl_Sust_all(~isnan(tbl_Sust_all.("COE [USD2010/kWh]")),:);
disp("Created long table w HP data")
%% Only box plot
figure;
subplot(2,2,2)
boxplot(PnetAllsComp,GrpscenptypComp,...%'orientation','horizontal',...
     'color',cmap3energy,'PlotStyle','compact','OutlierSize',2,'Symbol','.','FactorDirection','list');
grid on
set(gca,'YScale', 'log')
ylim([0.03, 10000])
addHPcategoryLines(1,'h')
title("PNet")
%
subplot(2,2,3)
boxplot(tbl_Sust{:,"COE [USD2010/kWh]"},tbl_Sust{:,"Scenario-plant"},'orientation','horizontal',...
     'color',cmap3energy,'PlotStyle','compact','OutlierSize',2,'FactorDirection','list'); %'label',{'',''},
xline(costlim,'LineStyle','--','Color','b','LineWidth',1.25,'DisplayName','Financial potential limit'); %Economic potential
grid on
set(gca,'XScale', 'log','XTickLabel','')
yline([2.5:2:2.5*2],'Color', mygraylines)
title("COE")

%% Only box plot - multipe grouping vars
figure;
boxplot(PnetAllsComp,{ScensComp SysIDsComp  },'color',cmap3energy,'whisker',inf,...'ColorGroup',ScensComp,...%'orientation','horizontal',... 
     'PlotStyle','compact','OutlierSize',2,'Symbol','.','FactorDirection','list','FullFactors','off','FactorGap',[30 ,.1],'FactorSeparator',[1]);
grid on
set(gca,'YScale', 'log')
ylim([0.03, 10000])
title("PNet")

%% GOOD: Manual scatter plot - grouped by Large, Medium, Small
grps=categories(tbl_Sust.("Scenario-plant"));
 figure
subplot(4,4,[1:3 5:7 9:11])
 hold all
 for i=1:2:length(grps)
    selrows=tbl_Sust.("Scenario-plant")==grps(i);
    s1=scatter(tbl_Sust{selrows,2},tbl_Sust{selrows,3},10+tbl_Sust{selrows,3}/8,cmap3energy(i,:),'^','filled','MarkerFaceAlpha',myalpha, 'MarkerEdgeAlpha',myalpha+.2, 'DisplayName',grps{i});
        selrows=tbl_Sust.("Scenario-plant")==grps(i+1);
    s2=scatter(tbl_Sust{selrows,2},tbl_Sust{selrows,3},10+tbl_Sust{selrows,3}/8,cmap3energy(i+1,:),'filled','MarkerFaceAlpha',myalpha,'MarkerEdgeAlpha',myalpha+.2, 'DisplayName',grps{i+1});
 end
    set(gca,'XScale', 'log','YScale', 'log','YAxisLocation','Right')

% Add annotations
ylim([0.01, 10000])
xline(costlim,'LineStyle','--','Color','b','LineWidth',1.25,'DisplayName','Financial potential limit'); %Economic potential
%text(costlim*1.15, 0, sprintf("\\it Economic potential\n <= %0.2f USD/kWh",costlim),'FontSize',12)
addHPcategoryLines(0,'h')
ylabel('GWh','FontWeight','bold')
xlabel('Per unit production cost (USD2010 per KWh)','FontWeight','bold')

l=legend([grps; 'Financial limit' ;'Hydropower Size Classes']  ); 
%l.Title.String="SCENARIO             PLANT TYPE";
applymyplotformat("3 energy Sust FULL case")
mainx=xlim();
mainy=ylim();

% Add right box plot
subplot(4,4,[4 8 12])
boxplot(tbl_Sust{:,"PnetAlls [GWh/yr]"},tbl_Sust{:,"Scenario-plant"},'label',{'','','','','','',},...%'orientation','horizontal',...
     'color',cmap3energy,'PlotStyle','compact','OutlierSize',2,'Symbol','.','FactorDirection','list');
ylim(mainy)
addHPcategoryLines(1,'h')
xline([2.5:2:2.5*2],'Color', mygraylines)
set(gca,'YScale', 'log','XTicklabel',[])
grid on

% Add bottom box plot
subplot(4,4,[13:15])
boxplot(tbl_Sust{:,"COE [USD2010/kWh]"},tbl_Sust{:,"Scenario-plant"},'orientation','horizontal',...
     'color',cmap3energy,'PlotStyle','compact','OutlierSize',2,'FactorDirection','list','label',{'','','','','','',})
xline(costlim,'LineStyle','--','Color','b','LineWidth',1.25,'DisplayName','Financial potential limit'); %Economic potential
xlim(mainx)
yline([2.5:2:2.5*2],'Color', mygraylines)
set(gca,'XScale', 'log','YTicklabel',[])
grid on

%% GOOD: Manual scatter plot - grouped by DP RP
grps=categories(tbl_Sust.("Scenario-plant"));
figure
subplot(4,4,[1:3 5:7 9:11])
 hold all
 
 for i=1:2:length(grps)
    selrows=tbl_Sust.("Scenario-plant")==grps(i);
    s1=scatter(tbl_Sust{selrows,2},tbl_Sust{selrows,3},10+tbl_Sust{selrows,3}/8,cmap3energy(i,:),'^','filled','MarkerFaceAlpha',myalpha, 'MarkerEdgeAlpha',myalpha+.2, 'DisplayName',grps{i});
        selrows=tbl_Sust.("Scenario-plant")==grps(i+1);
    s2=scatter(tbl_Sust{selrows,2},tbl_Sust{selrows,3},10+tbl_Sust{selrows,3}/8,cmap3energy(i+1,:),'filled','MarkerFaceAlpha',myalpha,'MarkerEdgeAlpha',myalpha+.2, 'DisplayName',grps{i+1});
 end
 
 
 set(gca,'XScale', 'log','YScale', 'log','YAxisLocation','Right')

% Add annotations
ylim([0.01, 10000])
xline(costlim,'LineStyle','--','Color','b','LineWidth',1.25,'DisplayName','Financial potential limit'); %Economic potential
%text(costlim*1.15, 0, sprintf("\\it Economic potential\n <= %0.2f USD/kWh",costlim),'FontSize',12)
addHPcategoryLines(0,'h')
ylabel('Energy (GWh per year)','FontWeight','bold')
xlabel('Per unit production cost (USD2010 per KWh)','FontWeight','bold')

l=legend([grps; 'Financial limit' ;'Hydropower size classes']  ); 
%l.Title.String="SCENARIO             PLANT TYPE";
applymyplotformat("3 energy Sust FULL case")
mainx=xlim();
mainy=ylim();

% Add right box plot
tbl_Sust{:,"Scenario-plant_2"}=categorical(tbl_Sust{:,"Scenario-plant"}, grps([2 4 6 1 3 5]));

subplot(4,4,[4 8 12])
boxplot(tbl_Sust{:,"PnetAlls [GWh/yr]"},tbl_Sust{:,"Scenario-plant_2"},... %'label',{'','','','','','',},'orientation','horizontal',...
     'color',cmap3energy([2 4 6 1 3 5],:),'PlotStyle','compact','OutlierSize',2,'Symbol','.','FactorDirection','list');
ylim(mainy)
addHPcategoryLines(1,'h')
xline(3.5,'Color', mygraylines)
set(gca,'YScale', 'log')%,'XTicklabel',[])
set(gca,'XTick',[1:6],'XTicklabel',extractBefore(grps([2 4 6 1 3 5]),"|"))
    grid on

% Add bottom box plot
subplot(4,4,[13:15])
boxplot(tbl_Sust{:,"COE [USD2010/kWh]"},tbl_Sust{:,"Scenario-plant_2"},'orientation','horizontal','label',{'','DP','','','RP','',},...  ... %
     'color',cmap3energy([2 4 6 1 3 5],:),'PlotStyle','compact','OutlierSize',2,'Symbol','.','FactorDirection','list');      
xline(costlim,'LineStyle','--','Color','b','LineWidth',1.25,'DisplayName','Financial potential limit'); %Economic potential
xlim(mainx)
yline(3.5,'Color', mygraylines)
set(gca,'XScale', 'log','XAxisLocation','Top')
set(gca,'YAxisLocation','Right','YTick',[1:6],'YTicklabel',extractBefore(grps([2 4 6 1 3 5]),"|"))
grid on


%% Split box plots RP and DP- not happy w it
figure;
selrows=ismember(tbl_Sust{:,"Scenario-plant_2"},grps([2,4,6]));
boxplot(tbl_Sust{selrows,"COE [USD2010/kWh]"},tbl_Sust{selrows,"Scenario-plant_2"},'orientation','horizontal',... %'label',{'','','','','','',},...
     'color',cmap3energy([2 4 6 1 3 5],:),'PlotStyle','compact','OutlierSize',5,'Symbol','o','FactorDirection','list');      
hold all

selrows=ismember(tbl_Sust{:,"Scenario-plant_2"},grps([1,3,5]));
boxplot(tbl_Sust{selrows,"COE [USD2010/kWh]"},tbl_Sust{selrows,"Scenario-plant_2"},'orientation','horizontal',... %'label',{'','','','','','',},...
     'color',cmap3energy([2 4 6 1 3 5],:),'PlotStyle','compact','OutlierSize',5,'Symbol','^','FactorDirection','list');      

xline(costlim,'LineStyle','--','Color','b','LineWidth',1.25,'DisplayName','Financial potential limit'); %Economic potential
xlim(mainx)
yline(3.5,'Color', mygraylines)
set(gca,'XScale', 'log','XAxisLocation','Top')
%set(gca,'YAxisLocation','Right','YTick',[1:6],'YTicklabel',grps([2 4 6 1 3 5]))
grid on


%% MEH: Box+Violin variations - Try gramm variations like in R
%g(1,1)=gramm('x',cars.Origin_Region,'y',cars.Horsepower,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5);
%g.set_names('x','Origin','y','Horsepower','color','# Cyl');

%Violin plots
g(1,1)=gramm('x',tbl_Sust.("Scen"),'y',tbl_Sust.("PnetAlls [GWh/yr]"),'color',tbl_Sust.("Scen"),'lightness',tbl_Sust.("Ptype"));
g(2,1)=copy(g(1));
g(2,2)=copy(g(1));

g(1,1).stat_violin('fill','transparent'); %'half',true,
g(1,1).set_title('stat_violin()');

%Boxplots+violin
g(1,2)=copy(g(1));
g(1,2).stat_boxplot();
g(1,2).set_title('stat_boxplot()');

% Jittered swarm
g(2,1).geom_jitter('width',0.6,'height',0,'dodge',0.7,'alpha',0.05);
g(2,1).set_title('jittered data');

%box+violin
 g(2,2).set_names('x','Search');
 g(2,2).stat_violin('dodge',0,'fill','edge');
 g(2,2).stat_boxplot('width',0.15);
% g(2,2).set_title('with stat_boxplot()');
% g(2,2).set_color_options('map','brewer_dark');

g(:).axe_property('YScale','log','YGrid','on' ) ;  
figure;g.draw();

%% MEH: Density histogram examples
% Kernel smoothing density estimate
clear g2
g2(1,1)=gramm('x',tbl_Sust.("PnetAlls [GWh/yr]"),'color',tbl_Sust.("Scen"),'lightness',tbl_Sust.("Ptype"));
g2(2,1)=copy(g2(1));

g2(1,1).facet_grid(tbl_Sust.("Ptype"),[]);
g2(1,1).stat_density();
g2(1,1).set_title('stat_density()');

% Transparent histogram
g2(2,1).facet_grid(tbl_Sust.("Ptype"),[]);
g2(2,1).stat_bin('fill','transparent');
g2(2,1).set_title('''transparent''');

g2(:).axe_property('XScale','log','XGrid','on' ) ;  
figure;g2.draw();
% didnt seem worth it to try this in matlab proper based on these trials but could follow this
% example: https://nl.mathworks.com/matlabcentral/answers/493279-greetings-anyway-to-help-to-create-ridgeline-by-matlab

%% 2D density contors - cant add this to a subplot - might be good to try this w the other params
% Display points and 95% percentile confidence ellipse
clear g2
g2(1,1)=gramm('x',tbl_Sust.("COE [USD2010/kWh]")    ,'y',tbl_Sust.("PnetAlls [GWh/yr]"),'color',tbl_Sust.("Scen"),'lightness',tbl_Sust.("Ptype"));

g2(1,1).geom_point('alpha',0.05);
g2(1,1).set_point_options('base_size',2);
g2(1,2)=copy(g2(1));
g2(1,1).set_names('x',"COE [USD2010/kWh]",'y',"PnetAlls [GWh/yr]",'color','Color=Search type','lightness','Lightness=Plant type');
g2(1,1).stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',2});
g2(1,1).set_title('stat_ellispe()');

%Plot point density as contour plot
g2(1,2).stat_bin2d('nbins',[10 10],'geom','contour');
g2(1,2).set_names('x',"COE [USD2010/kWh]",'y',"PnetAlls [GWh/yr]",'color','Color=Search type','lightness','Lightness=Plant type');
g2(1,2).set_title('stat_bin2d(''geom'',''contour'')');

g2(:).axe_property('YScale','log','XScale','log','XGrid','on' ) ;  
figure;g2.draw();
