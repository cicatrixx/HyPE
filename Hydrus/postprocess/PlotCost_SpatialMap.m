% Plot and create gif animating the cost curve for full mixed sust scenario
% Note that here we take the tech case of the sust scenarios.

clc
clearvars
close all
run('myVarNames.m')
load('G:\SurfDrive\HPmodel\output\Figs_trial\MainScenarios.mat', 'pcsout',...
    'basindata','catchments_cl')
animateGIF=0;

%% Sort data by cost from low to high
i=6; % the Full mixed sust scenario

COEAlls=pcsout{i}.COEAlls;
PnetAlls=pcsout{i}.PnetAlls;
coss=pcsout{i}.coss;
ross=pcsout{i}.ross;
SysIDs=pcsout{i}.SysIDs;

% Sort all cells in the order of the COEs from low to high
[COESort,ind]=sort(COEAlls);
% select only notnan inds
nanstart=find(isnan(COESort));
ind=ind(1:nanstart-1);
COESort=COESort(1:nanstart-1);

PnetSort=PnetAlls(ind); %in GWh
PnetSort_cum=cumsum(PnetSort)/1000; % convert GWh to TWh
%create symbol size vector
PnetSortSize= rescale(PnetSort,5,1100); %,'InputMax',maxscal,'InputMin',minscal);

SysIDsSort=SysIDs(ind);
coSort=coss(ind);
roSort=ross(ind);

% PnetISort_cum=[0; PnetISort_cum];
% COEISort=[0; COEISort];


%% Bubbleplot: Static Cost Curve vs Location of plant 
cmap8=cmap8_wong;
mycol=[1 1 1]*.3;

figure('Position', get(0, 'Screensize'))
subplot(3,3,[1,2,4,5])
h=imagescnan(catchments_cl);
set(h, 'AlphaData', myalpha*(catchments_cl>0))
colormap(cmap8)
hold on
bubblechart(coSort,roSort,PnetSortSize,mycol) %,"o","filled","MarkerFaceAlpha",0.3,"MarkerEdgeAlpha",0.9)
colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'Ticklabels',(basindata.basinnames),'TickDirection','none',...%  'Location','southoutside' ...
    'Position',    [0.175694444444444 0.363256784968684 0.402777777777778 0.0222686151704941])

subplot(3,3,[3,6])
scatter(PnetSort_cum,COESort,PnetSortSize,mycol) %,"o","MarkerFaceAlpha",myalpha)
hold on
stairs(PnetSort_cum,COESort,"-",'LineWidth',1.5,'Color',mycol,'DisplayName','Cost curve' );
xlabel('Cumulative Energy (TWh)')
ylabel('Per unit production cost (USD2010 per KWh)')
box on
grid on
set(gca,'YScale','log')

%%  Scatter: Static Cost Curve vs Location of plant w DP RP separate
i=6; % the Full mixed sust scenario
mycol=[1 1 1]*.3;

figure('Position', get(0, 'Screensize'))
subplot(1,6,1:4)
imagescnan(catchments_cl)
hold on
scatter(coSort(SysIDsSort==1),roSort(SysIDsSort==1),PnetSortSize(SysIDsSort==1),mycol,"o","filled","MarkerFaceAlpha",0.2,"MarkerEdgeColor",mycol);
scatter(coSort(SysIDsSort==2),roSort(SysIDsSort==2),PnetSortSize(SysIDsSort==2),mycol,"^","filled","MarkerFaceAlpha",0.2,"MarkerEdgeColor",mycol)
colormap(cmap8)
colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'Ticklabels',(basindata.basinnames),'TickDirection','none', 'Location','north')

subplot(1,6,5:6)
scatter(PnetSort_cum(SysIDsSort==1), COESort(SysIDsSort==1),PnetSortSize(SysIDsSort==1),mycol,"o","MarkerFaceAlpha",myalpha)
hold on
scatter(PnetSort_cum(SysIDsSort==2), COESort(SysIDsSort==2),PnetSortSize(SysIDsSort==2),mycol,"^","MarkerFaceAlpha",myalpha)
stairs(PnetSort_cum,COESort,"-",'LineWidth',1.5,'Color',mycol,'DisplayName','Cost curve' );
yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
text(300,costlim*1.15, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)
xlabel('Cumulative Energy (TWh)')
ylabel('Per unit production cost (USD2010 per KWh)')
box on
grid on
set(gca,'YScale','log','YTick',[0.01 0.1 1 10 50],'YTickLabel',[0.01 0.1 1 10 50])
set(gca,'FontName',    'Segoi UI Black','FontSize',14);

%% TRY Paused frames - OAT RP and DP
fig=figure('Position', get(0, 'Screensize'));
subplot(1,6,1:4)
imagescnan(catchments_cl)
colormap(cmap8)
hold on
bubblechart(coSort(1),roSort(1),PnetSortSize(1),mycol)

subplot(1,6,5:6)
stairs(PnetSort_cum,COESort,"-",'LineWidth',1.5,'Color',mycol,'DisplayName','Cost curve' );
hold on
yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
text(300,costlim*1.15, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)
xlabel('Cumulative Energy (TWh)')
ylabel('Per unit production cost (USD2010 per KWh)')
box on
grid on
set(gca,'YScale','log','YTick',[0.01 0.1 1 10 50],'YTickLabel',[0.01 0.1 1 10 50])

for idx=1:100%length(PnetSort_cum) % the Full mixed sust scenario
    if SysIDsSort(idx)==1
        subplot(1,6,1:4)
        scatter(coSort(idx),roSort(idx),PnetSortSize(idx),mycol,"o","filled","MarkerFaceAlpha",0.2,"MarkerEdgeColor",mycol)
        subplot(1,6,5:6)
        scatter(PnetSort_cum(idx), COESort(idx),PnetSortSize(idx),mycol,"o","filled","MarkerFaceAlpha",0.2,"MarkerEdgeColor",mycol)
    else
        subplot(1,6,1:4)
        scatter(coSort(idx),roSort(idx),PnetSortSize(idx),mycol,"^","filled","MarkerFaceAlpha",0.2,"MarkerEdgeColor",mycol)
        subplot(1,6,5:6)
        scatter(PnetSort_cum(idx), COESort(idx),PnetSortSize(idx),mycol,"^","filled","MarkerFaceAlpha",0.2,"MarkerEdgeColor",mycol)
    end
    sgtitle(idx)
    pause(0.01)
end


%% FINAL: Plot GIF of cost curve and map - no RP DP distinction
lims=[1:2:9 10:5:150 151:100:length(PnetSortSize) length(PnetSortSize)]; % looping through indices
nplots=length(lims);

fig=figure('Position', get(0, 'Screensize'),'color','w');
%subplot(1,6,1:4)
subtightplot(1,6,1:4, [0.02 0.06],[.12 .12],[.01 .01]); %[vgap hgap], marg_h,marg_w ;
h=imagescnan(catchments_cl);
set(h, 'AlphaData', myalpha*(catchments_cl>0))
colormap(cmap8)
hold on
%bubblechart(coSort(1),roSort(1),PnetSortSize(1),mycol)
cb=colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'Ticklabels',(basindata.basinnames),'TickDirection','none',  'Location','southoutside'); %
% Move colorbar lower so it doesnt change subplot size 
% position = [left bottom width height]
%cb.Position = cb.Position + 1e-10; 
cb.Position=[0.0099    0.1233    0.6333    0.02];

% Make the colorbar transparent - Need to do this manually
% cdata = cb.Face.Texture.CData;
% cdata(end,:) = uint8(myalpha * cdata(end,:));
% cb.Face.Texture.ColorType = 'truecoloralpha';
% cb.Face.Texture.CData = cdata;
set(gca,'FontName',    'Segoi UI Black','FontSize',14);

subtightplot(1,6,5:6, [0.02 0.06],[.16 .16],[.01 .02]); %[vgap hgap], marg_h,marg_w ;
stairs(PnetSort_cum,COESort,"-",'LineWidth',1.5,'Color',mycol,'DisplayName','Cost curve' );
hold on
yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
text(1,costlim*.85, sprintf("\\it Finanical potential <= %0.2f USD/kWh",costlim),'FontSize',12)
xlabel('Cumulative Energy (TWh)',    'FontWeight','bold')
ylabel('Per unit production cost (USD2010 per KWh)',    'FontWeight','bold')
box on
grid on
set(gca,'YScale','log','YTick',[0.01 0.1 1 10 50],'YTickLabel',[0.01 0.1 1 10 50])
set(gca,'FontName',    'Segoi UI Black','FontSize',14);

%%
selDPs=SysIDsSort==1;
selRPs=SysIDsSort==2;
for m=1:nplots-1
    %Save fig to frame
    drawnow
    frame = getframe(fig);
    im{m} = frame2im(frame);

    %if SysIDsSort(idx)==1
    idx=lims(m):lims(m+1);
    subtightplot(1,6,1:4, [0.02 0.06],[.12 .12],[.01 .01]); %[vgap hgap], marg_h,marg_w ;
    scatter(coSort(idx),roSort(idx),PnetSortSize(idx),mycol,"o","filled","MarkerFaceAlpha",0.2,"MarkerEdgeColor",mycol)%,"MarkerEdgeAlpha",0)
    subtightplot(1,6,5:6, [0.02 0.06],[.16 .16],[.01 .02]); %[vgap hgap], marg_h,marg_w ;
    scatter(PnetSort_cum(idx), COESort(idx),PnetSortSize(idx),mycol,"o","filled","MarkerFaceAlpha",0.2,"MarkerEdgeColor",mycol)%,"MarkerEdgeAlpha",0)
    %     else
    %         subplot(1,6,1:4)
    %         scatter(coSort(idx),roSort(idx),PnetSortSize(idx),mycol,"^","filled","MarkerFaceAlpha",0.2,"MarkerEdgeColor",mycol)
    %         subplot(1,6,5:6)
    %         scatter(PnetSort_cum(idx), COESort(idx),PnetSortSize(idx),mycol,"^","filled","MarkerFaceAlpha",0.2,"MarkerEdgeColor",mycol)
    %      end
    currenttitle=sprintf("Number of projects: %d\nTotal potential:  %0.0f TWh/yr",idx(end),PnetSort_cum(idx(end)));
% Create textbox
annotation(fig,'textbox',...
    [0.430729166666666 0.86 0.141666666666667 0.0605866388308843],...
    'String',currenttitle,...
    'LineStyle','none',...
    'FontSize',14,...
    'FontWeight','bold',...
    'FitBoxToText','off', ...
    'HorizontalAlignment','center', ...
    'BackgroundColor','w');


    pause(0.1)
    %Save fig to frame
    drawnow
    frame = getframe(fig);
    im{m} = frame2im(frame);
end
%Save fig to frame
drawnow
frame = getframe(fig);
im{m+1} = frame2im(frame);
%     %% Show plots as subplots
%     figure;
%     for m=1:nplots
%         subplot(ceil(nplots/2),2,m)
%         %figure
%         imshow(im{m});
%     end

%% Animate Plots
if animateGIF
    GIFpath=fullfile(rootf,'output','Figs_trial','cl');
    filename = fullfile(GIFpath,'GA_CostCurveVsSpatialMap.gif'); % Specify the output file name
    delaytsec=0.15;
    for selstations=1:nplots
        [A,map] = rgb2ind(im{selstations},256);
        if selstations == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delaytsec);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delaytsec);
        end
    end
end
