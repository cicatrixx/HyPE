%% Create plots comparing R and Q at diff outlets for observed and simulated
% downscaleQ_compare.m
% Created By    : Sanita Dhaubanjar on 02 Jun 2020
% Created For	: SustaIndus project WP2
% Last Updated  : 18 Mar 2020
%=========================

clc
close all
clear all
root=pwd; %"G:/GitConnect";
createplots=0;
wrapup=1;      % create potential figs for paper
implementkx=0; % implement SPHY routing
datapath=fullfile(root,"data","UI");
%matout=fullfile(fullfile(root,"output","downscaleQ","evalPerformance"), "RQatoutlets_m3s.mat");
figpath=fullfile(root,"output","TheoryPot_WrapUp");

%% Load 5km basinmask and outlets
data5km='UI5km.mat';
b5km=load(fullfile(datapath,data5km), 'outlet','catchment','outside');%,  'fdir', 'outlet'
b5km.outletID=unique(b5km.outlet(b5km.outlet~=0)); %~isnan(outlet)
countUniques(b5km.catchment)

% Load MAT5km Q
QMLAB="MLAB5km_40yrClimatology_m3day.mat";
load(fullfile(datapath,QMLAB),'','R_m3day');
b5km.Q_m3s_MLAB = double(cat(3,Q_m3day{:}))/(24*60*60);
b5km.Q_m3s_MLAB(b5km.outside)=nan;

% Load SPHY 5km Q and R maps
QSPHY = 'SPHY5km_40yrClimatology.mat';
load(fullfile(datapath,QSPHY));
b5km.Q_m3s_SPHY= double(cat(3,Qm_SPHY_m3sec{:}));
scaleunit=5000*5000*1e-3/(24*60*60); % convert mm/day to m3/s at 5km
%b5km.R_mmday = cat(3,Rm_SPHY_mmday{:});
b5km.R_m3s = double(cat(3,Rm_SPHY_mmday{:}))*scaleunit;

fprintf("Loaded 5km datasets \n\n")

%% Load downscaled 500m basinmask, catchments, outlets, Q and R maps
data500m='UI500m.mat';
b500m=load(fullfile(datapath,data500m), 'outlet','outlet_orig','catchment','outside_union','basinID','basinname','subbasinmask_union');%,'subbasinmask_matcat');%,  'fdir', 'outlet'
b500m.outside=b500m.outside_union;
b500m.outletID=unique(b500m.outlet_orig); %~isnan(outlet)
b500m.outletID=unique(b500m.outlet_orig(b500m.outlet_orig~=0)); %~isnan(outlet)
%
QMLAB="MLAB500m_40yrClimatology_m3day.mat";
load(fullfile(datapath,QMLAB),'Q500m_m3day','R500m_m3day');
b500m.Q_m3s_SPHY=cat(3,Q500m_m3day{:})/(24*60*60);
b500m.R_m3s=cat(3,R500m_m3day{:})/(24*60*60);

countUniques(b500m.catchment)
fprintf("Loaded 500m datasets \n\n")

all( b5km.outletID==b500m.outletID)
npts=length(b500m.outletID);

%% Mask basins to catchments
for m=1:13
    b5km.Q_m3s_SPHY(:,:,m)=(maskBasin(b5km.Q_m3s_SPHY(:,:,m),~b5km.outside));
    b500m.Q_m3s_SPHY(:,:,m)=(maskBasin(b500m.Q_m3s_SPHY(:,:,m),~b500m.outside));
end

%% FINAL: Plot 5km and 500m Q maps for 12 months
animateGIF=0;
if createplots
    %%
    nm=13;
    mname = month(datetime(2013,01,01):calmonths(1):datetime(2013,12,31),'name');
    mname = {mname{:},'ANNUAL'};
    Qmax =  max(max( b500m.Q_m3s_SPHY,[],'all'),max( b5km.Q_m3s_SPHY,[],'all'));
    %Set up fix for making GIF
    fig=figure('Position', get(0, 'Screensize'));
    axes('Units', 'normalized', 'Position', [0 0 1 1])
    for m=1:nm%:12
        subplot(1,2,1)
        imagescnan(log(b5km.Q_m3s_SPHY(:,:,m))); caxis([0, log(Qmax)])
        %imagescnan(b5km.Q_m3s_SPHY(:,:,m)); caxis([0, Qmax]); set(gca, 'ColorScale', 'log')
        axis image
        c=colorbar('Location','westoutside');
        c.Label.String = 'log(Discharge (m^3/s))';
        %c.Label.String = 'Discharge (m^3/s)';
        title('5km')
        
        subplot(1,2,2)
        imagescnan(log(b500m.Q_m3s_SPHY(:,:,m))); caxis([0, log(Qmax)])
        %imagescnan(b500m.Q_m3s_SPHY(:,:,m)); caxis([0, Qmax]);set(gca, 'ColorScale', 'log')
        axis image
        
        title('500m')
        sgtitle(mname{m})
        
        drawnow
        frame = getframe(fig);
        im{m} = frame2im(frame);
    end
    
    %% Show plots
    %     figure;
    %     for idx = 1:nm
    %         subplot(ceil(nm/2),2,idx)
    %         imshow(im{idx});
    %     end
    
    %% animate plots
    if animateGIF
        filename = fullfile(figpath,'Q5kmVs500m_logval.gif'); % Specify the output file name
        delaytsec=1.25;
        for selstations = 1:nm
            [A,map] = rgb2ind(im{selstations},256);
            if selstations == 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delaytsec);
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delaytsec);
            end
        end
    end
    %orient('landscape');print('tmp.pdf','-dpdf','-bestfit')
end

%% Compare Qs: Extract station TS from Q5km
mydata=b5km;
Q5km_TS=zeros(14,npts,'double');
for pt=1:length(mydata.outletID)
    % Get 5km TS
    [r0,c0]=find(mydata.outlet==mydata.outletID(pt));
    Q5km_TS(1:14,pt) = squeeze(mydata.Q_m3s_SPHY(r0,c0,:));
end

%% Compare Qs: Extract station TS from Q500m
mydata=b500m;
Q500m_TS=zeros(14,npts,'double');
for pt=1:length(mydata.outletID)
    % Get 500m TS
    [r0,c0]=find(mydata.outlet_orig==mydata.outletID(pt));
    Q500m_TS(1:13,pt) = squeeze(mydata.Q_m3s_SPHY(r0,c0,:));
end
fprintf("Extracted TS from 5km and 500m Qmaps\n\n")
%% Compare Qs: Create comparison for Q5km vs Q500m TS - all outlets
cmap = linspecer(npts); %
%cmap = cmap(randperm(size(cmap,1)), :); %shuffle color order
if createplots
    %% Plot catchments and outlets w IDs
    tmpcat=b500m.catchment;
    %only in catchment matrix change discontinuous indices for correct
    %colormapping
    newID=[8 10 12:15; 6:11];
    for i = 1:length(newID)
        tmpcat(tmpcat==newID(1,i)) = newID(2,i);
    end
    
    tmp=find(b500m.outlet);
    ptID=b500m.outlet(tmp);
    [ir0,sc]=ind2sub(size(b500m.outlet), tmp);
    
    %b500m.catchment(b500m.catchment==-1)=nan;
    figure
    imagescnan(tmpcat)
    text(sc-50, ir0+100, num2str(ptID),'FontSize' ,8,'BackgroundColor',.9*[1 1 1]);%,'FontWeight','bold'
    hold all; plot(sc,ir0,".k",'markersize',15)   %plot outlets
    colormap([0 0 0; cmap])
    colorbar('XTick', 0:1:11,'XTickLabel',num2str(b500m.outletID))
    grid on
    title("Outlet IDs and Catchments")
    
    %% Monthly TS plot
    figure;
    h1=plot(Q5km_TS(1:12,:),'LineWidth',1);
    xlabel("Months")
    ylabel("Q in m^3/s")
    xlim([1 12])
    grid on
    hold on
    
    %figure
    h2=plot(Q500m_TS(1:12,:),'*', 'LineWidth',2,'markersize',8);
    set(h1, {'color'}, num2cell(cmap,2));
    set(h2, {'color'}, num2cell(cmap,2));
    set(gca, 'YScale', 'log')
    legend( num2str(mydata.outletID(2:end)),'Location','EastOutside')
    title("Monthly Q simulated by SPHY(-) and Matlab(*) ")
    
    %% Create crossplots
    figure
    %Monthly
    subplot(1,2,1)
    line([0 8000],[0 8000],'Color',[.5 .5 .5],'LineWidth',1,'LineStyle',":")
    hold on
    for i=1:npts
        plot(Q5km_TS(1:12,i),Q500m_TS(1:12,i),"*",'Color',cmap(i,:),'MarkerSize',8);
    end
    xlabel("Q5km")
    ylabel("Q500m")
    axis square
    grid on
    title("Monthly averages at 11 points")
    legend(["x=y line";  num2str(mydata.outletID(2:end))],'Location','SouthEast')
    
    %Annual
    subplot(1,2,2)
    line([0 3000],[0 3000],'Color',[.5 .5 .5],'LineWidth',1,'LineStyle',":")
    hold on
    for i=1:npts
        plot(Q5km_TS(13,i),Q500m_TS(13,i),"*",'Color',cmap(i,:),'MarkerSize',10);
        plot(Q5km_TS(14,i),Q500m_TS(13,i),"o",'Color',cmap(i,:),'MarkerSize',10);
    end
    xlabel("Q5km")
    ylabel("Q500m")
    axis square
    grid on
    title("Annual averages at 11 points")
    
end

%% Compare w Runoffs: Extract station TS as sum R5km
mydata=b5km;
R5km_TS=zeros(14,npts,'double');
for pt=1:length(mydata.outletID)
    catcells=mydata.catchment==mydata.outletID(pt);
    % Get 5km sum(Runoff in cat) for each tstep
    for t=1:13
        R5km_TS(t,pt) = nansum(maskBasin(mydata.R_m3s(:,:,t),catcells),'all');
    end
end

% if all( R5km_TS(:)==Q5km_TS(:))
%     disp("At 5km, sum(R5km) do not match Q(5km)")_
% end

%% Compare w Runoffs:  Extract station TS as sum R500m
mydata=b500m;
R500m_TS=zeros(14,npts,'double');
for pt=1:length(mydata.outletID)
    catcells=mydata.catchment==mydata.outletID(pt);
    % Get 5km sum(Runoff in cat) for each tstep
    for t=1:13
        R500m_TS(t,pt) = nansum(maskBasin(mydata.R_m3s(:,:,t),catcells),'all');
    end
end

%% Eval perfstats for R and Qs
% Monthly
[monPerf.Pbias(1,:), monPerf.NSE(1,:)]= evalStats(R5km_TS(1:12,:), Q5km_TS(1:12,:));
[monPerf.Pbias(2,:), monPerf.NSE(2,:)]= evalStats(R5km_TS(1:12,:), Q500m_TS(1:12,:));
[monPerf.Pbias(3,:), monPerf.NSE(3,:)]= evalStats(R5km_TS(1:12,:), R500m_TS(1:12,:));
[monPerf.Pbias(4,:), monPerf.NSE(4,:)]= evalStats(R500m_TS(1:12,:), Q500m_TS(1:12,:));

if createplots
    %% Plot: Differences in R and Q
    figure
    subplot(2,2,1);imagesc( R5km_TS-Q5km_TS); title("sum(R5km)-Qkm");caxis([1e-6 1500])
    xlabel("Stations"); ylabel("12 months + annual")
    subplot(2,2,2);imagesc( R5km_TS-Q500m_TS); title("sum(R500m)-Q500m");caxis([1e-6 1500])
    xlabel("Stations"); ylabel("12 months + annual")
    subplot(2,2,3);imagesc( R5km_TS-R500m_TS); title("sum(R5km)-Qkm");caxis([1e-6 1500])
    xlabel("Stations"); ylabel("12 months + annual")
    subplot(2,2,4);imagesc( R500m_TS-Q500m_TS); title("sum(R500m)-Q500m");caxis([1e-6 1500]);
    xlabel("Stations"); ylabel("12 months + annual")
    c=colorbar;
    c.Label.String = "m^3/s";
    sgtitle("Differences between outlet-wise sumR and Q at 5km and 500m")
    
    %% Summary Crossplot: SPHY sum(R5km) vs SPHY Q5km, MLAB Q500m and sum(R500m)
    figure;
    plot(R5km_TS(1:12,:),Q5km_TS(1:12,:),"s")
    hold all
    plot(R5km_TS(1:12,:),Q500m_TS(1:12,:),"o")
    plot(R5km_TS(1:12,:),R500m_TS(1:12,:),".")
    lgd=legend([repmat(" ",22,1);num2str(mydata.outletID(2:end))],...
        'Location', 'BestOutside','NumColumns', 3);
    lgd.Title.String="Q5km        Q500m        R500m StnID";
    %set(gca, 'YScale', 'log','XScale', 'log')
    xlabel("Station-wise sum of R5km (m^3/s)")
    ylabel("Processed (m^3/s)")
    title("Comparison of SPHY R5km with SPHY Q5km and Matlab R500m and Q500m")
    %
    applymyplotformat(cmap)
    %grid minor
    line([0 8000],[0 8000],'Color',[.5 .5 .5],'LineWidth',1,'LineStyle',":")
end
fprintf("Compared performance of R and Q \n")

%% Compare w observed Q and eval perfstats
selpts = find(ismember(b500m.outletID,[3 8 10]));
obsQ = readmatrix(fullfile(datapath, "obsQ.xlsx"),'Sheet','Obs_climatology','Range','B3:D14');

[monPerf2.Pbias(1,:), monPerf2.NSE(1,:)]= evalStats(obsQ(1:12,:), Q5km_TS(1:12,selpts));
[monPerf2.Pbias(2,:), monPerf2.NSE(2,:)]= evalStats(obsQ(1:12,:), Q500m_TS(1:12,selpts));
[monPerf2.Pbias(3,:), monPerf2.NSE(3,:)]= evalStats(obsQ(1:12,:), R5km_TS(1:12,selpts));
[monPerf2.Pbias(4,:), monPerf2.NSE(4,:)]= evalStats(obsQ(1:12,:), R500m_TS(1:12,selpts));

%% FINAL: Catchment and stations for perfEval
cmap = linspecer(length(b500m.basinID),'qualitative');
selcs=find(ismember(b500m.basinID,[103 104 105]));
selstations=find(ismember(b500m.outlet,[103 104 105]));
stationnames=["Tarbela", "Mangala", "Marala"];
[sr,sc]=ind2sub(size(b500m.outlet), selstations);
% Plot catchments w checked station IDs
idx=find(~isnan(b500m.outlet) & b500m.outlet ~=110); %find only outlets that matter
ptID=b500m.outlet(idx);
[ir0,ic0]=ind2sub(size(b500m.outlet), idx);

f=figure('Position', get(0, 'Screensize'));
%Plot all subbasin and labels
imagescnan(b500m.subbasinmask_union)
text(ic0, ir0,categorical(ptID,b500m.basinID,b500m.basinname),'FontSize' ,6)%,'BackgroundColor',.9*[1 1 1]);%,'FontWeight','bold'
hold all;
%plot selected stations and station IDs
plot(sc,sr,".k",'markersize',15)   
text(sc, sr,stationnames,'FontSize' ,8,'FontWeight','bold') %'BackgroundColor',.9*[1 1 1]
colormap(cmap)
c=colorbar('Ticks', b500m.basinID,'TickLabels',b500m.basinname,'Direction','reverse');
c.Label.String="Sub-basins";
grid on
title("Outlet IDs and Catchments")
%orient('landscape'); print('G:\GitConnect\output\TheoryPot_WrapUp\UnionCatchments_QcheckStations.pdf','-dpdf','-fillpage')

%% FINAL: Plot obs vs sim Q for select stations only
figure;
% Monthly Cross plot
subplot(1,2,2)
hold all
plot(obsQ(1:12,:),Q5km_TS(1:12,selpts),"+",'markersize',6, 'LineWidth',1.2)
plot(obsQ(1:12,:),Q500m_TS(1:12,selpts),"x",'markersize',8, 'LineWidth',1.2)
line([0 8000],[0 8000],'Color','k','LineWidth',1,'LineStyle',":")
xlabel("Observed (m^3/s)")
ylabel("Simulated (m^3/s)")
applymyplotformat(cmap(selcs,:))
%grid minor
lgd=legend([repmat(" ",3,1); "   Tarbela"; "   Mangala" ;"   Marala"],...
    'Location', 'NorthOutside','NumColumns', 2);
lgd.Title.String="Sim-5km  Sim-500m   StnID";


% Monthly TS plot
subplot(1,2,1)
plot(obsQ(1:12,:),'.-','LineWidth',1.2);
hold all
plot(Q5km_TS(1:12,selpts),'+:', 'LineWidth',1.2,'markersize',6);
plot(Q500m_TS(1:12,selpts),'x:', 'LineWidth',1.2,'markersize',8);
lgd=legend([repmat(" ",3*2,1); "   Tarbela"; "   Mangala" ;"   Marala"],...
    'Location', 'NorthOutside','NumColumns', 3);
lgd.Title.String="   Obs           Sim-5km     Sim-500m    StnID  ";
applymyplotformat(cmap(selcs,:))
xlabel("Months")
ylabel("Q in m^3/s")
xlim([1 12])

sgtitle("Monthly observed and simulated Q/R in SPHY and Matlab")
%orient('landscape'); print('G:\GitConnect\output\TheoryPot_WrapUp\Qsimvsobs.pdf','-dpdf','-fillpage')

%% Plot sim vs observed - everything

if createplots
    %% Summary Crossplot:
    figure;
    subplot(1,2,2)
    hold all
    plot(obsQ(1:12,:),Q5km_TS(1:12,selpts),"+",'markersize',6)
    plot(obsQ(1:12,:),Q500m_TS(1:12,selpts),"s",'markersize',8)
    plot(obsQ(1:12,:),R5km_TS(1:12,selpts),"o",'markersize',12)
    plot(obsQ(1:12,:),R500m_TS(1:12,selpts),".",'markersize',10)
    line([0 8000],[0 8000],'Color','k','LineWidth',1,'LineStyle',":")
    xlabel("Observed (m^3/s)")
    ylabel("Simulated (m^3/s)")
    %
    applymyplotformat(linspecer(3))
    %grid minor
    
    %% Monthly TS plot
    subplot(1,2,1)
    plot(obsQ(1:12,:),'LineWidth',1.2);
    hold all
    plot(Q5km_TS(1:12,selpts),'+:', 'LineWidth',1,'markersize',6);
    plot(Q500m_TS(1:12,selpts),'s:', 'LineWidth',1,'markersize',8);
    plot(R5km_TS(1:12,selpts),'o', 'LineWidth',1,'markersize',12);
    plot(R500m_TS(1:12,selpts),'.:', 'LineWidth',1,'markersize',10);
    lgd=legend([repmat(" ",3*4,1); "Tarbela"; "Mangala" ;"Marala"],...
        'Location', 'NorthOutside','NumColumns', 5);
    lgd.Title.String="Obs Q5km-SPHY    Q500m-MLAB       R5km-SPHY    R500m-MLAB   StnID";
    
    applymyplotformat(linspecer(3))
    xlabel("Months")
    ylabel("Q in m^3/s")
    xlim([1 12])
    sgtitle("Monthly observed and simulated Q/R in SPHY and Matlab")
    fprintf("Compared performance of simulated and observed \n")
    
end

%% Create individual crossplots of R and Q --old
if createplots
    %% Crossplot: SPHY sum(Rtot5km) vs SPHY Q5km
    figure
    %Monthly
    subplot(1,2,1)
    line([0 8000],[0 8000],'Color',[.5 .5 .5],'LineWidth',1,'LineStyle',":")
    hold on
    for i=1:npts
        plot(R5km_TS(1:12,i),Q5km_TS(1:12,i),"*",'Color',cmap(i,:),'MarkerSize',8);
    end
    xlabel("SPHY sum R5km")
    ylabel("SPHY Q5km")
    axis square
    grid on
    title("Monthly averages at 11 points")
    legend(["x=y line";  num2str(mydata.outletID(2:end))],'Location','SouthEast')
    
    %Annual
    subplot(1,2,2)
    line([0 3000],[0 3000],'Color',[.5 .5 .5],'LineWidth',1,'LineStyle',":")
    hold on
    for i=1:npts
        plot(R5km_TS(13,i),Q5km_TS(13,i),"*",'Color',cmap(i,:),'MarkerSize',10);
        %plot(R5km_TS(14,i),Q500m_TS(13,i),"o",'Color',cmap(i,:),'MarkerSize',10);
    end
    xlabel("SPHY sum R5km")
    ylabel("SPHY Q5km")
    axis square
    grid on
    title("Annual averages at 11 points")
    sgtitle("SPHY sum R5km vs SPHY Q 5km")
    
    %% Crossplot: SPHY sum(Rtot5km) vs Matlab Q500m
    figure
    %Monthly
    subplot(1,2,1)
    line([0 8000],[0 8000],'Color',[.5 .5 .5],'LineWidth',1,'LineStyle',":")
    hold on
    for i=1:npts
        plot(R5km_TS(1:12,i),Q500m_TS(1:12,i),"*",'Color',cmap(i,:),'MarkerSize',8);
    end
    xlabel("SPHY sum R5km")
    ylabel("MATLAB Q500m")
    axis square
    grid on
    title("Monthly averages at 11 points")
    legend(["x=y line";  num2str(mydata.outletID(2:end))],'Location','SouthEast')
    
    %Annual
    subplot(1,2,2)
    line([0 3000],[0 3000],'Color',[.5 .5 .5],'LineWidth',1,'LineStyle',":")
    hold on
    for i=1:npts
        plot(R5km_TS(13,i),Q500m_TS(13,i),"*",'Color',cmap(i,:),'MarkerSize',10);
        %plot(R5km_TS(14,i),Q500m_TS(13,i),"o",'Color',cmap(i,:),'MarkerSize',10);
    end
    xlabel("SPHY sum R5km")
    ylabel("MATLAB Q500m")
    axis square
    grid on
    title("Annual averages at 11 points")
    sgtitle("SPHY sum R5km vs MATLAB Q500m")
    %% Crossplot: Matlab sum(Rtot500m) vs Matlab  Q500m
    figure
    %Monthly
    subplot(1,2,1)
    line([0 8000],[0 8000],'Color',[.5 .5 .5],'LineWidth',1,'LineStyle',":")
    hold on
    for i=1:npts
        plot(R500m_TS(1:12,i),Q500m_TS(1:12,i),"*",'Color',cmap(i,:),'MarkerSize',8);
    end
    xlabel("MATLAB sum R500m")
    ylabel("MATLAB Q500m")
    axis square
    grid on
    title("Monthly averages at 11 points")
    legend(["x=y line";  num2str(mydata.outletID(2:end))],'Location','SouthEast')
    
    %Annual
    subplot(1,2,2)
    line([0 3000],[0 3000],'Color',[.5 .5 .5],'LineWidth',1,'LineStyle',":")
    hold on
    for i=1:npts
        plot(R500m_TS(13,i),Q500m_TS(13,i),"*",'Color',cmap(i,:),'MarkerSize',10);
        %plot(R5km_TS(14,i),Q500m_TS(13,i),"o",'Color',cmap(i,:),'MarkerSize',10);
    end
    xlabel("MATLAB sum R500m")
    ylabel("MATLAB Q500m")
    axis square
    grid on
    title("Annual averages at 11 points")
    sgtitle("SPHY sum R5km vs MATLAB Q500m")
end

%% Implement kx shifting in 5km
if implementkx
    kx=0.95;        %from Sonu
    Qaccu=b5km.Q_m3s_SPHY;
    Qrout=Qaccu*0;
    
    for t=1:12
        if t==1
            Qrout(:,:,t)=Qaccu(:,:,t);%(1-kx)*Qaccu(:,:,t)+ kx*Qrout; % for jan take leftover flow from dec
        else
            Qrout(:,:,t)=(1-kx)*Qaccu(:,:,t)+ kx*Qrout(:,:,t-1);
        end
    end
    b5km.Q_m3s_rout=Qrout;
    fprintf("Shifted Q5km \n")
    
    %% Plot Q maps
    m=12;
    res="5km";
    figure
    subplot(1,3,1)
    imagesc(b5km.Q_m3s_SPHY(:,:,m))
    set(gca, 'ColorScale', 'log')
    colorbar
    axis image
    title('SPHY')
    
    subplot(1,3,2)
    imagesc(b5km.Q_m3s_MLAB(:,:,m))
    set(gca, 'ColorScale', 'log')
    colorbar
    axis image
    title('Gernaat')
    
    subplot(1,3,3)
    imagesc(Qrout(:,:,m))
    set(gca, 'ColorScale', 'log')
    colorbar
    axis image
    title('Gernaat + kx')
    sgtitle(sprintf("%s at m=%d",res,m))
end
%% save
if exist('matout')
    save(matout, 'Q5km_TS' ,'Q500m_TS', 'R5km_TS' ,'R500m_TS','MonPerf')
end
disp("***************************************EOF***************************************")

