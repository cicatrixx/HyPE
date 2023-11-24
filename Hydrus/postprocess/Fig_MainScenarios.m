%% Make barplots for the 3 main energy scenarios, the 3 geohazards and the 3 geohazard risk scenarios.
clc
clearvars
close all
addpath(genpath('G:\SurfDrive\GitConnect\Hydrus\devFiles\'), ...
    genpath('G:\SurfDrive\HPmodel\Hydrus\'))
run('myVarNames.m')

load(fullfile(pwd,'output','Compile_Scenario_Totals03.mat'))
cl=tot_allscen{:,1};
figure
%t=tiledlayout(2,4,'Padding','compact','TileSpacing','compact');

ylimitsallplots=[0 320];

%% Add empty rows to force space between three energy scens
tot_allscen{96,1}={'x'};
tot_allscen{96,2:7}=nan;
%% FINAL: Plot bar for 3 energy scenarios full and remain overlaid
%nexttile([1 4])
subplot(2,4,[1:4])
selfull=[1:3 96 4:6 96 7:9] ; %1:9; %mainidx_end;
selremain=[10:12 96 13:15 96 16:18]; %mainidx_end;
b1=bar(tot_allscen{selfull, 2:3},.95,'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',baralpha,'EdgeAlpha',0);
hold all
b2=bar(tot_allscen{selremain, 2:3},.7,  'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'LineStyle','-');
ylim(ylimitsallplots)

% add lines separating boxes
%    xline([3.5,6.5],'Color', 'k')
ylabel('Energy in TWh/year','FontWeight','bold')
i=1;text(0.01,300,char(64+i) ,'FontWeight','bold','FontSize',16)


% Change color of  bar
for k = 1:2
    b1(k).FaceColor = mycolors{k};
    b1(k).EdgeColor = mycolors{k};
    
    b2(k).EdgeColor = mycolors{k}*0.6;
end
yyaxis right
scshift=0.1;
scatter([1:length(selfull)]-scshift,tot_allscen{selfull, 5},30,cl_RP*.6,'d','filled')
hold on
scatter([1:length(selfull)]-scshift,tot_allscen{selfull, 6}+tot_allscen{selfull, 5},30, cl_DP*.6,'d','filled') %manually stack scatter
%
scatter([1:length(selfull)]+scshift,tot_allscen{selremain, 5},30,cl_RP*.6,'d','LineWidth',1.15)
scatter([1:length(selfull)]+scshift,tot_allscen{selremain, 6}+tot_allscen{selremain, 5},30, cl_DP*.6,'d','LineWidth',1.15) %manually stack scatter
% fix yyaxis
%ylabel('\diamondsuit = Number of projects','FontWeight','bold')
ylabel('Number of projects','FontWeight','bold')
ax = gca;
ax.YAxis(2).Color = 'k';

% add scenario labels
xticklabels(repmat([pottypes3 {' '}],1,3))
xlabel([""; "Energy Focus Scenarios"],'FontWeight','bold')
text([2 5 8],-800*[1 1 1],searchtypes,'HorizontalAlignment','center','FontWeight','bold','FontAngle','italic')
% add legend
l=legend(["Full", "Full", "Remain |", "Remain |","Full","Full", "Remain", "Remain"], "NumColumns", 4,'Location','northoutside'); %,'Orientation','horizontal')%)
l.Title.String=sprintf("Energy Potential | Number of Projects ");

% legend(["RP: Full Energy", "DP: Full Energy", "RP: Remaining Energy", "DP: Remaining Energy",...
%     "RP: Full Number","DP: Full Number", "RP: Remaining Number","DP: Remaining Number"], "NumColumns", 4,'Location','northoutside') %,'Orientation','horizontal')%)
%box off; 
 grid on
set(gca,'XGrid','off')
% Create textbox w RP and DP labels
annotation('textbox',...
    [0.118589743129472 0.906335452019717 0.183194149339996 0.0590094822224375], ...    %[0.184516444185945 0.852821668148887 0.215419501133787 0.0999999999999999],...
    'VerticalAlignment','bottom',...
    'String',['River power plant:',newline(),'Diversion canal plant:'],...
    'HorizontalAlignment','right',...
    'FontSize',10,...
    'FontName','segoe ui',...
    'FontWeight','bold',...
    'FitBoxToText','on',...    
    'LineStyle','none');
%% FINAL: Plot waterfall bar
%nexttile([1 2])
subplot(2,4,5:6);
selidx=[7 8 25 24 9]; %25
plotWaterFallDecline2Data(tot_allscen{selidx, 2}', tot_allscen{selidx, 3}', tech2sust_Labels, [cl_RP; cl_DP])
ylim(ylimitsallplots)
xtickangle(40)
xlabel('Succesive addition of constraints','FontWeight','bold')
ylabel('Energy in TWh/year','FontWeight','bold')
set(gca,'XGrid','off')
%box off; 
i=i+1;text(0.01,300,char(64+i) ,'FontWeight','bold','FontSize',16)

%% FINAL: Plot bar for geo hazard type and geo hazard rep scenarios
%figure
%nexttile
subplot(2,4,7)
selidx=[19:21,9];
b1=bar(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 2:3},'stacked','FaceColor',"flat",'FaceAlpha',baralpha,'EdgeAlpha',0);
grid on; set(gca,'XGrid','off')
ylabel('Energy in TWh/year','FontWeight','bold')
ylim(ylimitsallplots)
i=i+1;text(0.01,300,char(64+i) ,'FontWeight','bold','FontSize',16)

yyaxis right
scatter(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 5},30,cl_RP*.6,'d','filled')
hold on
scatter(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 6}+tot_allscen{selidx, 5},30, cl_DP*.6,'d','filled') %manually stack scatter
%Setup
%legend(t.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
%ylabel('\diamondsuit = Number of projects','FontWeight','bold')
xticklabels([geohazard_names "All 3"]); xtickangle(30)
xlabel("Geo-hazard types",'FontWeight','bold')
% Change color of  bar
for k = 1:2
    b1(k).FaceColor = mycolors{k};
    b1(k).EdgeColor = mycolors{k};
end
alpha(b1,0.7)
ax = gca;
ax.YAxis(2).Color = 'k';
set(gca,'YTick',[])
%box off; 

%nexttile
subplot(2,4,8)
selidx=[9,22:23];
b1=bar(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 2:3},'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',baralpha,'EdgeAlpha',0);
set(gca,'YTickLabels',"")
ylim(ylimitsallplots)

%ylabel('Energy in TWh/year','FontWeight','bold')
yyaxis right
scatter(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 5},30,cl_RP*.6,'d','filled')
hold on
scatter(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 6}+tot_allscen{selidx, 5},30, cl_DP*.6,'d','filled') %manually stack scatter

%Setup
%legend(t.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
ylabel('\diamondsuit = Number of projects','FontWeight','bold')
xticklabels(geohazard_scennames_cl([3,2,4])); xtickangle(40)
xlabel(['Geo-hazard risk', newline(), 'representations'],'FontWeight','bold')
% Change color of  bar
for k = 1:2
    b1(k).FaceColor = mycolors{k};
    b1(k).EdgeColor = mycolors{k};
end
alpha(b1,0.7)
ax = gca;
ax.YAxis(2).Color = 'k';
grid on
set(gca,'XGrid','off')
%box off; 
i=i+1;text(0.01,300,char(64+i) ,'FontWeight','bold','FontSize',16)

%export_fig(fullfile(pwd,'output','Figs_trial','cl','F3_MainScenarios'))