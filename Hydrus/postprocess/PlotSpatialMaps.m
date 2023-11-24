% Load data for main scenarios compiled in PlotSubBasinWiseData code
% - same as gettotals
clc
clearvars
close all
addpath(genpath('G:\SurfDrive\GitConnect\Hydrus\devFiles\'), ...
    genpath('G:\SurfDrive\HPmodel\Hydrus\'))

run('myVarNames.m')

load('G:\SurfDrive\HPmodel\output\Figs_trial\MainScenarios4R.mat', 'pcsout',...
    'runnames','pcsouthaz',  'runnameshaz','basindata','catchments_cl')
 pcsout=[pcsout pcsouthaz];
runnames=[runnames runnameshaz];

ifname = sprintf('%s\\data\\%s\\Basin_UIB\\PantpeBasin_%d.mat', rootf, continent_in,nbasin);
datain=load(ifname,'outside','channel_main_trib');
cmap8=brighten(cmap8_wong,0);

%% GOOD: Loop to make 2 col spatial maps w optimal sites under Tech/Econ/Sust 3 energy scenarios
runnames_cl=erase(strrep(extractAfter(runnames,'Full_'),'_','-'),["-Fin"]);
myscens=[1:6];
nplot=length(myscens);

%get basin outline
inbasin=single(~datain.outside);
inbasin(datain.outside)=nan;

% get channel
[rch,cch]=ind2sub(size(datain.outside), find(datain.channel_main_trib));
minmksz=5;

%create symbol size vector
%mydata.PnetAllSize= rescale(PnetSort,25,1500); %,'InputMax',maxscal,'InputMin',minscal);

figure
nplot=length(myscens);
fig_order= 1:nplot;      % [5,6,3,4,1,2];%   [1,3,5,2,4,6]; %
for i= 1:length(myscens) %
    mydata=pcsout{myscens(i)};
    if myscens(i)==1; mydata.co_arc=[];mydata.ro_arc=[];end %older run where arc was not saved
%    subplot(nplot/2,2,fig_order(i))
                subtightplot (nplot/2,2,fig_order(i), [0.025 0.025]); %[vgap hgap] subplot(3,3,tt);

    plottitle=runnames_cl{myscens(i)}; %""; %
    plotMyRaster(inbasin,'',plottitle);axis image;colormap([1 1 1]*0.9);colorbar off

    % plot subbasins
%     h=imagescnan(catchments_cl);hold all
%     colormap(cmap8)
%     set(h, 'AlphaData', myalpha*(catchments_cl>0))
%     title(plottitle)

hold on
    l1=scatter(cch,rch,.5,cl_channel,'.','DisplayName',"Channel"); %plot channel
    %l2=plot([mydata.co_arc; mydata.cos],[mydata.ro_arc; mydata.ros],'.','Color',mycolor_reject,'markersize',mksz-2,'DisplayName',"Rejected Sites"); %,'Color',0.7*[1 1 1]
    l4=scatter(mydata.coss(mydata.SysIDs==2),mydata.ross(mydata.SysIDs==2),sqrt(mydata.PnetAlls(mydata.SysIDs==2))+minmksz,cl_RP,'^','filled','MarkerEdgeAlpha',0.4,'DisplayName',planttypes{1}); 
    l3=scatter(mydata.coss(mydata.SysIDs==1),mydata.ross(mydata.SysIDs==1),sqrt(mydata.PnetAlls(mydata.SysIDs==1))+minmksz,cl_DP*.6,'o','MarkerEdgeAlpha',0.4,'DisplayName',planttypes{2}); 
    text(100,100,char(64+i),'FontName','franklin gothic medium','FontSize',14)   % add label ABC to top-left of fig
    %set(gca,'FontName',    'franklin gothic medium','FontSize',14,'Layer','top');
end
legend([l1 l4 l3 ],'location','bestoutside','Orientation','horizontal')
%legend([l3 l4 l2 l1],'location','bestoutside','Orientation','horizontal')

% savefig(fullfile(outfig_fldr,'OptimalLocs_AllMajor.fig'))
%savefig(fullfile(outfig_fldr,'OptimalLocs2_ForEGU2021.fig'))
%export_fig(fullfile(outfig_fldr,'OptimalLocs2_ForEGU2021.jpg'))

%% GOOD: Loop to make 2 col spatial maps w optimal sites under Tech/Econ/Sust 3 HAZARD REP SCENARIOS
runnames_cl=erase(strrep(extractAfter(runnames,'Full_'),'_','-'),["Mixed-"])';
myscens=[6 13 14];
nplot=length(myscens);

%get basin outline
inbasin=single(~datain.outside);
inbasin(datain.outside)=nan;

% get channel
[rch,cch]=ind2sub(size(datain.outside), find(datain.channel_main_trib));
minmksz=5;

%create symbol size vector
%mydata.PnetAllSize= rescale(PnetSort,25,1500); %,'InputMax',maxscal,'InputMin',minscal);

figure
nplot=length(myscens);
fig_order= 1:nplot;      % [5,6,3,4,1,2];%   [1,3,5,2,4,6]; %
for i= 1:length(myscens) %
    mydata=pcsout{myscens(i)};
    if myscens(i)==1; mydata.co_arc=[];mydata.ro_arc=[];end %older run where arc was not saved
%    subplot(nplot/2,2,fig_order(i))
                subtightplot (2,2,fig_order(i), [0.025 0.025]); %[vgap hgap] subplot(3,3,tt);

    plottitle=runnames_cl{myscens(i)}; %""; %
    plotMyRaster(inbasin,'',plottitle);axis image;colormap([1 1 1]*0.9);colorbar off
    hold on
    l1=scatter(cch,rch,1,cl_channel,'.','DisplayName',"Channel"); %plot channel
    %l2=plot([mydata.co_arc; mydata.cos],[mydata.ro_arc; mydata.ros],'.','Color',mycolor_reject,'markersize',mksz-2,'DisplayName',"Rejected Sites"); %,'Color',0.7*[1 1 1]
    l4=scatter(mydata.coss(mydata.SysIDs==2),mydata.ross(mydata.SysIDs==2),sqrt(mydata.PnetAlls(mydata.SysIDs==2))+minmksz,cl_RP,'^','filled','MarkerEdgeAlpha',0.4,'DisplayName',planttypes{1}); 
    l3=scatter(mydata.coss(mydata.SysIDs==1),mydata.ross(mydata.SysIDs==1),sqrt(mydata.PnetAlls(mydata.SysIDs==1))+minmksz,cl_DP*.6,'o','MarkerEdgeAlpha',0.4,'DisplayName',planttypes{2}); 
    text(100,100,char(64+i),'FontName','franklin gothic medium','FontSize',14)   % add label ABC to top-left of fig
    %set(gca,'FontName',    'franklin gothic medium','FontSize',14,'Layer','top');
end
legend([l1 l4 l3 ],'location','bestoutside','Orientation','horizontal')
%legend([l3 l4 l2 l1],'location','bestoutside','Orientation','horizontal')

% savefig(fullfile(outfig_fldr,'OptimalLocs_AllMajor.fig'))
%savefig(fullfile(outfig_fldr,'OptimalLocs2_ForEGU2021.fig'))
%export_fig(fullfile(outfig_fldr,'OptimalLocs2_ForEGU2021.jpg'))
