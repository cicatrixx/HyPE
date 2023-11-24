%% Create plots comparing Q at diff outlets for observed and simulated (Q vs Qwc)
% downscaleQ_compare.m
% Created By    : Sanita Dhaubanjar on 7 Dec 2020
% Created For	: SustaIndus project WP2
%=========================

% Check if Qwc created correctly
clc
close all
clear all
root=pwd; %"G:/GitConnect";
compareQnatQwc=1;
compareQobsQsimNat=0;
wrapup=1;      % create potential figs for paper
implementkx=0; % implement SPHY routing
datapath=fullfile(root,"data","UI","data");
figpath=fullfile(root,"output","downscaleQwc2");
matout=fullfile(figpath, "RQatoutlets_m3s_R4.mat");

%% Load 5km basinmask and outlets
data5km='UI5km.mat';
b5km=load(fullfile(datapath,data5km), 'outlet','catchment','outside');%,  'fdir', 'outlet'
b5km.outletID=unique(b5km.outlet(b5km.outlet~=0)); %~isnan(outlet)
% disp("Catchment cells in 5km")
% countUniques(b5km.catchment)

% Load MAT5km Q
QMLAB="MLAB5km_40yrClimatology_m3day.mat";
load(fullfile(datapath,QMLAB),'Q_m3day','R_m3day');
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

%% Load downscaled 500m basinmask, catchments, outlet_orig, Q and R maps
data500m='UI500m_ArcGIS.mat';

b500m=load(fullfile(datapath,data500m), 'outlet_orig','catchments','outside','basinlabels', 'outlets');%,'subbasinmask_matcat');%,  'fdir'
b500m.catchment=b500m.catchments;
b500m.basinID=b500m.basinlabels.basinIDs;
b500m.basinname=b500m.basinlabels.basinnames;
b500m.outlet=b500m.outlets;

% b500m=load(fullfile(datapath,data500m), 'outlet','outlet_orig','catchment','outside_union','basinID','basinname','subbasinmask_union');%,'subbasinmask_matcat');%,  'fdir', 'outlet'
% b500m.outside=b500m.outside_union;
b500m.outletID=unique(b500m.outlet_orig(b500m.outlet_orig~=0)); %~isnan(outlet)
%
QMLAB="MLAB500m_40yrClimatology_m3day_R4.mat";
load(fullfile(datapath,QMLAB),'Q500m_m3day','R500m_m3day','Qwc_500m_m3day','Rwc_500m_m3day');
b500m.Q_m3s_SPHY=cat(3,Q500m_m3day{:})/(24*60*60);
b500m.R_m3s=cat(3,R500m_m3day{:})/(24*60*60);
b500m.Qwc_m3s_SPHY=cat(3,Qwc_500m_m3day{:})/(24*60*60);
b500m.Rwc_m3s=cat(3,Rwc_500m_m3day{:})/(24*60*60);
fprintf("Loaded 500m datasets \n\n")

%% Mask basins
all( b5km.outletID==b500m.outletID(1:11))
npts=length(b500m.outletID);

for m=1:13
    b5km.Q_m3s_SPHY(:,:,m)=(maskBasin(b5km.Q_m3s_SPHY(:,:,m),~b5km.outside));
    b500m.Q_m3s_SPHY(:,:,m)=(maskBasin(b500m.Q_m3s_SPHY(:,:,m),~b500m.outside));
    b500m.Qwc_m3s_SPHY(:,:,m)=(maskBasin(b500m.Qwc_m3s_SPHY(:,:,m),~b500m.outside));
end


%% FINAL: Plot 5km and 500m Q and Qwc maps for 12 months
animateGIF=0;
if compareQnatQwc
    %%
    nm=13;
    mname = month(datetime(2013,01,01):calmonths(1):datetime(2013,12,31),'name');
    mname = {mname{:},'ANNUAL'};
    Qmax =  max(max( b500m.Q_m3s_SPHY,[],'all'),max( b5km.Q_m3s_SPHY,[],'all'));
    %Set up fix for making GIF
    fig=figure('Position', get(0, 'Screensize'));
    axes('Units', 'normalized', 'Position', [0 0 1 1])
    for m=1:nm%:12
        sgtitle(mname{m})
        
        subplot(2,2,1:2)
        imagescnan(log(b5km.Q_m3s_SPHY(:,:,m))); caxis([0, log(Qmax)])
        %imagescnan(b5km.Q_m3s_SPHY(:,:,m)); caxis([0, Qmax]); set(gca, 'ColorScale', 'log')
        axis image
        c=colorbar('Location','westoutside');
        c.Label.String = 'log(Discharge (m^3/s))';
        %c.Label.String = 'Discharge (m^3/s)';
        title('5km: SPHY Q in m^3/s')
        
        subplot(2,2,3)
        imagescnan(log(b500m.Q_m3s_SPHY(:,:,m))); caxis([0, log(Qmax)])
        %imagescnan(b500m.Q_m3s_SPHY(:,:,m)); caxis([0, Qmax]);set(gca, 'ColorScale', 'log')
        axis image
        title('500m: Matlab Q in m^3/s')
        
        subplot(2,2,4)
        imagescnan(log(b500m.Qwc_m3s_SPHY(:,:,m))); caxis([0, log(Qmax)])
        %imagescnan(b500m.Qwc_m3s_SPHY(:,:,m)); caxis([0, Qmax]);set(gca, 'ColorScale', 'log')
        axis image
        
        title('500m: Matlab Qwc in m^3/s')
        
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
        filename = fullfile(figpath,'Q5kmVsQ500m_Qwc500m_logval.gif'); % Specify the output file name
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
    disp('Animation of R, Q and Qwc saved')
end

%% Get TS from Qobs, Qnat and Qwc at all old outlets
for getQ_TS=1
    %% Compare Qs: Extract station TS from Q5km - uses outletIDs because new basin IDs not assigned to 5k
    mydata=b5km;
    Q5km_TS=zeros(14,npts,'double');
    for pt=1:length(mydata.outletID)
        % Get 5km TS
        [r0,c0]=find(mydata.outlet==mydata.outletID(pt));
        Q5km_TS(1:14,pt) = squeeze(mydata.Q_m3s_SPHY(r0,c0,:));
    end
    
    %% Compare Qs: Extract station Q TS from Q500m- uses outletIDs because new basin IDs not assigned to 5km
    mydata=b500m;
    Q500m_TS=zeros(14,npts,'double');
    for pt=1:length(mydata.outletID)
        % Get 500m TS
        [r0,c0]=find(mydata.outlet_orig==mydata.outletID(pt));
        Q500m_TS(1:13,pt) = squeeze(mydata.Q_m3s_SPHY(r0,c0,:));
    end
    
    %% Compare Qs: Extract station Qwc TS from Q500m- uses outletIDs because new basin IDs not assigned to 5km
    mydata=b500m;
    Qwc_500m_TS=zeros(14,npts,'double');
    for pt=1:length(mydata.outletID)
        % Get 500m TS
        [r0,c0]=find(mydata.outlet_orig==mydata.outletID(pt));
        Qwc_500m_TS(1:13,pt) = squeeze(mydata.Qwc_m3s_SPHY(r0,c0,:));
    end
    fprintf("Extracted TS from 5km and 500m Qmaps\n\n")
end

%% Compare Qnat vs Qwc
if compareQnatQwc
    %% FINAL: Plot monthly TS for Q500m vs Qwc500m  - all outlets
    oldbasinnames=b500m.basinname([2 1 3:end]); %Swat and kabul are switched from old to new for better figure ordering
    subplot = @(m,n,p) subtightplot (m, n, p, [0.003 0.02]);
    
    validpts_idx=[1,2,3,6,7,9,10,11]; %Exclude col numbers for the 110 equivalent old ids
    validbasinnames=string(categorical(b500m.outletID(validpts_idx), b500m.basinlabels.old_basinIDs,b500m.basinlabels.basinnames));
    nvpts=length(validpts_idx);
    cmap = linspecer(nvpts); %3 are 110 id
    
    figure;
    subplot(1,2,1)
    h1=plot(Q500m_TS(1:12,validpts_idx),'LineWidth',1);
    xlabel("Months")
    ylabel("Q in m^3/s")
    xlim([1 12])
    grid on
    hold on
    
    %figure
    h2=plot(Qwc_500m_TS(1:12,validpts_idx),'*', 'LineWidth',2,'markersize',8);
    set(h1, {'color'}, num2cell(cmap,2));
    set(h2, {'color'}, num2cell(cmap,2));
    set(gca, 'YScale', 'log')
    lgd=  legend([repmat(" ",nvpts,1);validbasinnames],'NumColumns', 2);  %num2str(ptID(validpts))
    lgd.Title.String="Qnat     Qwc    StnID";
    title("Monthly Qnat(line) and Qwc(*) simulated by Matlab ")
    
    % Monthly % reduction in Q
    subplot(1,2,2)
    h1=plot((Q500m_TS(1:12,validpts_idx)-Qwc_500m_TS(1:12,validpts_idx))./Q500m_TS(1:12,validpts_idx)*100,'*:');%,'LineWidth',1);
    xlabel("Months")
    ylabel("% reduction in Q after water consumption")
    legend(b500m.basinname),
    set(h1, {'color'}, num2cell(cmap,2));    xlim([1 12])
    grid on
    hold on
    title('Difference between Q and Qwc')
    
    % savefig(fullfile(figpath,'Q500m_Qwc500m_monthlyTS.fig'))
    
    %% FINAL: compare Qnatural and Qwc for annual
    basinmask=~b500m.outside;
    m=13;
    Qdiff=(Q500m_m3day{m}-Qwc_500m_m3day{m})./Q500m_m3day{m}*100;
    
    figure;
    subplot(2,2,1)
    imagescnan(maskBasin(Q500m_m3day{m},basinmask))
    mycbar('Qnatural in m^3/day')
    set(gca, 'ColorScale', 'log')
    
    subplot(2,2,2)
    imagescnan(maskBasin(Qwc_500m_m3day{m},basinmask))
    mycbar('Qwc in m^3/day')
    set(gca, 'ColorScale', 'log')
    
    subplot(2,2,3:4)
    imagescnan(maskBasin(Qdiff,basinmask))%'(Qnatural-Qwc)/Qnatural*100')
    mycbar('Percentage reduction in Q after WC','NorthOutside')
    
    sgtitle("For Annual")
    % savefig(fullfile(figpath,'QnatvsQwc_annualspatial.fig'))
    
end

%% Compare w observed and simulated Q  and eval perfstats for 3 stns
if compareQobsQsimNat
    selpts = find(ismember(b500m.outletID,[3 8 10]));
    obsQ = readmatrix(fullfile(datapath, "obsQ.xlsx"),'Sheet','Obs_climatology','Range','B3:D14');
    [Pbias(1,:), NSE(1,:)]= evalStats(obsQ(1:12,:), Q5km_TS(1:12,selpts));
    [Pbias(2,:), NSE(2,:)]= evalStats(obsQ(1:12,:), Q500m_TS(1:12,selpts));
    [Pbias(3,:), NSE(3,:)]= evalStats(obsQ(1:12,:), Qwc_500m_TS(1:12,selpts));
    Qobs_vs={'Qnat-5km'; 'Qnat-500m'; 'Qwc-500m'};
    monPerf=table(Qobs_vs, Pbias, NSE);
    % Catchment and selected stations for perfEval
    cmap = linspecer(length(b500m.basinID),'qualitative');
    selcmaps=find(ismember(b500m.basinID,[103 104 105]));
    stationnames=["Tarbela", "Mangala", "Marala"];
    
    %% Plot loc of outlets for sel stations
    selstations=find(ismember(b500m.outlet,[103 104 105]));
    [sr,sc]=ind2sub(size(b500m.outlet), selstations);
    
    
    %Plot all catchments w labels for selected station IDs
    %f=figure('Position', get(0, 'Screensize'));
    figure
    %Plot all subbasin
    imagescnan(maskBasin(b500m.catchments,~b500m.outside))
    hold all;
    %plot selected stations and station IDs
    plot(sc,sr,".k",'markersize',15)
    text(sc*1.01, sr,stationnames,'FontSize' ,18,'FontWeight','bold') %'BackgroundColor',.9*[1 1 1]
    colormap(cmap)
    c=colorbar('Ticks', b500m.basinID,'TickLabels',b500m.basinname,'Direction','reverse');
    c.Label.String="Sub-basins";
    grid on
    title("Outlet IDs and Catchments")
    
    %savefig(fullfile(figpath,'LocOfQVerification.fig'))
    
    %% FINAL: Plot obs vs sim Qnat for 3 selected stations only
    figure;
    % Monthly Cross plot
    subplot(1,2,2)
    hold all
    plot(obsQ(1:12,:),Q5km_TS(1:12,selpts),"+",'markersize',6, 'LineWidth',1.2)
    plot(obsQ(1:12,:),Q500m_TS(1:12,selpts),"x",'markersize',8, 'LineWidth',1.2)
    line([0 8000],[0 8000],'Color','k','LineWidth',1,'LineStyle',"--")
    xlabel("Observed (m^3/s)")
    ylabel("Simulated (m^3/s)")
    applymyplotformat('Qobs vs Qsim',cmap(selcmaps,:))
    
    %grid minor
    %legend(["   Tarbela"; "   Mangala" ;"   Marala"])
    
    % if Q and Qwc
    %plot(obsQ(1:12,:),Qwc_500m_TS(1:12,selpts),"+",'markersize',8, 'LineWidth',1.2)
    % lgd=legend([repmat(" ",3,1); "   Tarbela"; "   Mangala" ;"   Marala"],...
    %     'Location', 'NorthOutside','NumColumns', 2);
    % lgd.Title.String="Qnat-Sim     Qwc-Sim    StnID";
    
    % Monthly TS plot
    subplot(1,2,1)
    plot(obsQ(1:12,:),'.-','LineWidth',1.2);
    hold all
    plot(Q5km_TS(1:12,selpts),'+--', 'LineWidth',1.2,'markersize',6);
    plot(Q500m_TS(1:12,selpts),'x:', 'LineWidth',1.2,'markersize',8);
    applymyplotformat('Qobs vs Qsim',cmap(selcmaps,:))
    xlabel("Months")
    ylabel("Q in m^3/s")
    xlim([1 12])
    
    lgd=legend([repmat(" ",3*2,1); "   Tarbela"; "   Mangala" ;"   Marala"],...
        'Location', 'NorthOutside','NumColumns', 3);
    lgd.Title.String="  Qobs           Q5km     Q500m    StnID  ";
    
    
    % if Q and Qwc
    %plot(Qwc_500m_TS(1:12,selpts),'+-.', 'LineWidth',1.2,'markersize',8);
    %lgd=legend([repmat(" ",3*2,1); "   Tarbela"; "   Mangala" ;"   Marala"],...
    %    'Location', 'NorthOutside','NumColumns', 3);
    %lgd.Title.String="  Qobs           Qnat-Sim     Qwc-Sim    StnID  ";
    
    sgtitle("Monthly observed and simulated Q/Qwc in SPHY and Matlab")
    
    %orient('landscape'); print('G:\GitConnect\output\TheoryPot_WrapUp\Qsimvsobs.pdf','-dpdf','-fillpage')
    % savefig(fullfile(figpath,'QobsVSq5kmQ500m_3stations.fig'))
    
end

%% Compare annual Q with 15s data and validation points
load('G:\PaperData\DavidGernaat\Sanita_model_package\data\ASIA\Basin\Basin_5.mat', 'Q')
% ArcGIS based catchments
nodataval=-2147483648; %32bitsigned
ifldr='500m';
validpt_UI = loadSPHYtiff(fullfile(pwd,'data','data_prep',ifldr,"valid_pts_UIprj.tif"), nodataval , 0);
validpt_15sI = loadSPHYtiff(fullfile(pwd,'data','data_prep',ifldr,"valid_pts_David.tif"), nodataval , 0);
fprintf("Loaded 15s David dataset and valid pts \n\n")

%% Compare annual avg Qs: Extract at 3 valid pts
validID=[3,8,10];
nVpts=length(validID);
for pt=1:3
    % Get 15s val
    [r0,c0]=find(validpt_15sI==validID(pt));
    Q15s_TS(1,pt) = Q(r0,c0);
    % Get 500m val
    [r1,c1]=find(validpt_UI==validID(pt));
    Q500m_TS(1,pt) = squeeze(b500m.Qwc_m3s_SPHY(r1,c1,13));
end
 mean(obsQ)
%% plot difference in annual means for 3 station in the indus
figure;bar([mean(obsQ); Q500m_TS; Q15s_TS]')
xticklabels(["Tarbela", "Mangala" , "Marala"])
xlabel("Discharge stations")
ylabel("Historical long term annual average in m^3/s")
l=legend('Observed','Simulated: Khanal et al. 2020','Simulated: Gernaat et al. 2017');
grid on
%% save
if exist('matout')
    save(matout, 'Q5km_TS' ,'Q500m_TS','monPerf')%'R5km_TS' ,'R500m_TS'
    disp('Saved TS and perfstats')    
end
disp("***************************************EOF***************************************")
    
