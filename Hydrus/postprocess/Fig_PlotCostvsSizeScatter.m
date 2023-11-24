% Plot multi-panel graph w scatter and box plots. Scatter plot has symbol
% size scaled by plant size
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


%% Create color for 3 energy scenarios
%cmap3_mainenergy=cbrewer2('Set2',3);  % dark2 or set2 are color blind and print friendly
cmap3=[];
for i=1:3 % 3 scens % create lighter shade for each scenario
    % RP is lighter and plotted first
    cmap3=[cmap3; cmap3_mainenergy(i,:) ; brighten(cmap3_mainenergy(i,:),-0.6)];
end
% reorder and lighten colors a bit for box plots
cmap3_lighter=brighten(cmap3([2,4,6,1,3,5],:),0.4); % for box plots as they dont have opacity
%% Create table of 3 energy Sust 
myscens=[2,4,6];  const="Sust";%Sust case
%myscens=[1,3,5]; const="Tech"; %Tech case - need to change axis limits for

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
tbl_Sust_all = table(categorical(GrpscenptypComp,grps),COEAllsComp,PnetAllsComp,categorical(SysIDsComp,flip(planttypes)),categorical(ScensComp,searchtypes),...
    'VariableNames',{ 'Scenario-plant','COE [USD2010/kWh]', 'PnetAlls [GWh/yr]','Ptype','Scen'});
%select only not nans
tbl_Sust=tbl_Sust_all(~isnan(tbl_Sust_all.("COE [USD2010/kWh]")),:);
disp("Created long table w HP data")

%% FINAL: Manual scatter plot - grouped by DP RP
loglabels=[0.01, 0.1 1 10 100 1000 10000 100000];
grps=categories(tbl_Sust.("Scenario-plant"));
if const=="Sust"
    mainy=[0.009 10000];
else
    mainy=[0.009 100000];
end
mainx=[0.009 1000];

f=figure;
s(1)=subplot(4,4,[1:3 5:7 9:11]);
hold all
mymkrline=[0,0,2,2,0.5,0.5];
for i=1:2:length(grps)
    selrows=tbl_Sust.("Scenario-plant")==grps(i);
     if i>1
        s1=scatter(tbl_Sust{selrows,2},tbl_Sust{selrows,3},10+tbl_Sust{selrows,3}/8,cmap3(i,:),'^','MarkerFaceAlpha',myalpha, 'MarkerEdgeAlpha',myalpha+.1, 'DisplayName',grps{i},'LineWidth',mymkrline(i));
        selrows=tbl_Sust.("Scenario-plant")==grps(i+1);
        s2=scatter(tbl_Sust{selrows,2},tbl_Sust{selrows,3},10+tbl_Sust{selrows,3}/8,cmap3(i+1,:),'o','MarkerFaceAlpha',myalpha,'MarkerEdgeAlpha',myalpha+.1, 'DisplayName',grps{i+1},'LineWidth',mymkrline(i));
     elseif i==1
        s1=scatter(tbl_Sust{selrows,2},tbl_Sust{selrows,3},10+tbl_Sust{selrows,3}/8,cmap3(i,:),'filled','^','MarkerFaceAlpha',myalpha, 'MarkerEdgeAlpha',0, 'DisplayName',grps{i});
        selrows=tbl_Sust.("Scenario-plant")==grps(i+1);
        s2=scatter(tbl_Sust{selrows,2},tbl_Sust{selrows,3},10+tbl_Sust{selrows,3}/8,cmap3(i+1,:),'filled','o','MarkerFaceAlpha',myalpha,'MarkerEdgeAlpha',0, 'DisplayName',grps{i+1});
    end
end

set(gca,'XScale', 'log','YScale', 'log')  
set(gca,'YAxisLocation','Right')
xline(costlim,'LineStyle','--','Color','b','LineWidth',1.25,'DisplayName','Financial potential limit'); %Economic potential
addHPcategoryLines(01,'h',50)

%Fix log numbers on X and Y axis
yticks(loglabels)
xticks(loglabels)
yticklabels(num2str(yticks'))
xticklabels(num2str(xticks')) %use number labels w full numbers and not the exponential form!
ylim(mainy); xlim(mainx)

% Add annotations
ylabel('Energy (GWh per year)','FontWeight','bold')
xlabel('Per unit production cost (USD2010 per KWh)','FontWeight','bold')
%text(costlim*1.15, 0, sprintf("\\it Economic potential\n <= %0.2f USD/kWh",costlim),'FontSize',12)
%l=legend([grps; 'Financial limit' ;'Hydropower size classes'],'Location','southeastoutside'  );   %,'NumColumns',2
l=legend([" "; " ";" "; " "; planttypes{1};planttypes{2};"Financial limit" ;"Hydropower size classes"],'Location','northoutside' ,'NumColumns',7,'Box','off');
title(l,sprintf(" %s   %s    %s        PLANT TYPE               ",searchtypes{:}))
t=title(sprintf("3 energy %s FULL case",const));
grid on 
box on

% Add right box plot
s(2)=subplot(4,4,[4 8 12]);
boxplot(tbl_Sust{:,"PnetAlls [GWh/yr]"},tbl_Sust{:,["Ptype", "Scen"] },'color',cmap3_lighter, 'whisker',inf,'label',{},  ...'ColorGroup',ScensComp,...%'orientation','horizontal',...
    'PlotStyle','compact','OutlierSize',2,'Symbol','.','FactorDirection','data','FullFactors','off','FactorSeparator',[1]); %'FactorGap',[10 ,.1],
grid on
addHPcategoryLines(0,'h')
set(gca,'YScale', 'log','XGrid','off','FontSize',10) %,'XTicklabel',[])
%ylabel('Energy (GWh per year)','FontWeight','bold')
planttypes_2{2}=["River";"Power"; "Plant"];
planttypes_2{1}=["Diversion";"Canal"; "Plant"];
text([2,5],.02*[1 1],planttypes_2,'HorizontalAlignment','center','FontSize',10)
ylim(mainy)
set(gca,'YTicklabel',[],'XTick',[1:6],'XTicklabel',extractBefore(grps([2 4 6 1 3 5]),"|"))

% Add bottom box plot
s(3)=subplot(4,4,[13:15]);
boxplot(tbl_Sust{:,"COE [USD2010/kWh]"},tbl_Sust{:,["Ptype", "Scen"] },'color',cmap3_lighter, 'whisker',inf,'label',{},'orientation','horizontal',  ...'ColorGroup',ScensComp,...%'orientation','horizontal',...
    'PlotStyle','compact','OutlierSize',2,'Symbol','.','FactorDirection','list','FullFactors','off','FactorSeparator',[1]); %'FactorGap',[10 ,.1],
grid on
xline(costlim,'LineStyle','--','Color','b','LineWidth',1.25,'DisplayName','Financial potential limit'); %Economic potential
set(gca,'XScale', 'log','YGrid','off','YAxisLocation','Right','FontSize',10)%,'YTicklabel',[],'XAxisLocation','Top')
%xlabel('Per unit production cost (USD2010 per KWh)','FontWeight','bold')
text(500*[1 1],[2,5],planttypes_2,'HorizontalAlignment','center','FontSize',10)
xlim(mainx)
set(gca,'XTicklabel',[],'YTick',[1:6],'YTicklabel',extractBefore(grps([2 4 6 1 3 5]),"|"))


% Align the axes positions
set(s(1),'Position', [ 0.1          0.3          0.6          0.6]);
set(s(2),'Position', [ 0.76         0.3          0.12         0.6]);
set(s(3),'Position', [ 0.1          0.13          0.6          0.12]);
set(t,'Position',[30.64      70000             0]);
set(l,'Position',[0.12          0.91          0.57          0.07]);

annotation(f,'textbox',...
    [0.80 0.15 0.0814196242171218 0.0769768177028448],...
    'String',{'Energy Focus Scenarios'},...
    'FontWeight','bold',...
     'LineStyle','none');

