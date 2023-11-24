close all
clear all
%% Load data and select scenario
load('G:\SurfDrive\HPmodel\output\Figs_trial\MainScenarios_BasinOutAll.mat')

%% Load output data
%selrun= 'Test_Full_Mixed_Tech_Econ';
scen='Full';    % 'Remain';%
continent_in = 'ASIA';%RDVolumeLake15s
nbasin=103;
selrun= 'R103_Energy_Full_Mixed_Sust_RiskAverse'; 
root_outdata = fullfile(pwd,'output',scen,continent_in,selrun);
ofname=strcat(root_outdata,'\\Basin103_output_do.mat');

load(ofname)
disp("Loaded output data")
%
TarbelaVol_MCM=13.69*1000^3/1e6; % from wiki, max dam reservoir in the world is lake volta 8502km2 or 8502000MCM
TarbelaHeight_m=143.26; % from wiki, max in the world is 305m

cmap8=cmap8_wong;
%% Load input data
ifname = sprintf('%s\\data\\%s\\Basin_UIB\\PantpeBasin_%d.mat', root, continent_in, nbasin);
load(ifname,'acc', 'fdir','Z','Q','AllWaterBodies','outside','Regions');

run('defdirs.m')
catchments_cl=maskBasin(Regions,~outside);

%% Get data for valid RD points
PID  = find(SysID==1);
RDID = find(SysID==2);

ros_rdid = ross(RDID);
cos_rdid = coss(RDID);

valid_rdid=find(~isnan(ros_rdid)); % This is same for cos_rdid

rrd = ros_rdid(valid_rdid);
crd = cos_rdid(valid_rdid);
RDDepths = RDDepth(valid_rdid);

OptDHs_m = OptDH(valid_rdid); %in m
RDPnets_GWh = RDPnet(valid_rdid); % in kWh in functions but converted to GWh in main script
ResVols_MCM = RDVolumeLake15s(valid_rdid)/10e6; % m3 convert to MCM
Qdesigns_m3s=  Q_RD_design(valid_rdid);
COEs= COETotRD(valid_rdid);

% Sort RD projects by power size
[~, idx]=sort(RDPnets_GWh,'descend');

% Sort RD projects by res size
[~, idx_MCM]=sort(ResVols_MCM,'descend');

%% Load or scan reservoir surfaces
if loadresdata
    load(fullfile(root_outdata,'Reservoirs'),'RD_allSum','RDlake_Opt','RDlake','Q_burnR')
else
%% Trace the reservoirs - time consuming
parfor i=1:length(idx)
    fprintf('%d van %d\n',i,length(idx))
    RDinlet(i) = sub2ind(size(acc),rrd(idx(i)),crd(idx(i)));
    RDupstream = fastfindupstream_lim(acc,fdir,drow,dcol,RDinlet(i),0);
    RDlake_Opt{i}= i*(RDupstream & Z < (single(Z(RDinlet(i)))-RDDepths(idx(i))+OptDHs_m(idx(i))));
end

% plot just the lakes - time consuming
RDlake=cat(3,RDlake_Opt{:});
RD_allSum=sum(RDlake,3);
RD_allSum(RD_allSum==0)=nan;

figure;
imagescnan(RD_allSum)
hold on

%% Create Q overlaid w new reservoirs
Q_burnR= truecolorsc(log(Q),flipud(gray));
for i=1:numel(RDlake_Opt)
    Q_burnR = burnmask(Q_burnR, ~RDlake_Opt{i});
end

%% Save run as it is timeconsuming
save(fullfile(root_outdata,'Reservoirs'),'RD_allSum','RDlake_Opt','RDlake','Q_burnR','-v7.3')
end

%% Plot subbasins+newres
figure
h0=imagescnan(catchments_cl);
set(h0, 'AlphaData', .5*(catchments_cl>0))
hold on
tmp1=108+(RD_allSum>0);
tmp1(isnan(RD_allSum))=nan;
h1=imagescnan(tmp1); % all nan vals are set to transparent in this function
colormap([ cmap8; mlabblue;])
cb=colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'TickLabels',basinlabels.basinnames(selbasins),'TickDirection','none','Limits',[101 108] );%'Location','southoutside',
cb=colorbar('Ticks', [101.4 102.2 103 103.8 105.4 106.1 107 107.7],'TickLabels',basindata.basinnames,'TickDirection','none','Limits',[101 108.2] );%'Location','southoutside',

%     % Make the colorbar transparent - Need to do this manually
%     cdata = cb.Face.Texture.CData;
%     cdata(end,:) = uint8(myalpha * cdata(end,:));
%     cb.Face.Texture.ColorType = 'truecoloralpha';
%     cb.Face.Texture.CData = cdata;


%% Plot subbasins+newres+ reservoirs by symbolsize
figure
h0=imagescnan(catchments_cl);
set(h0, 'AlphaData', .5*(catchments_cl>0))
hold on
tmp1=108+(RD_allSum>0);
tmp1(isnan(RD_allSum))=nan;
h1=imagescnan(tmp1); % all nan vals are set to transparent in this function
colormap([ cmap8; mlabblue;])
cb=colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'TickLabels',basinlabels.basinnames(selbasins),'TickDirection','none','Limits',[101 108] );%'Location','southoutside',
%     % Make the colorbar transparent - Need to do this manually
%     cdata = cb.Face.Texture.CData;
%     cdata(end,:) = uint8(myalpha * cdata(end,:));
%     cb.Face.Texture.ColorType = 'truecoloralpha';
%     cb.Face.Texture.CData = cdata;

bubblechart(crd,rrd,ResVols_MCM,'k' ,'DisplayName',"Reservoir size indicators",'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.4) %,'MarkerEdgeColor',[1 1 1]*.2)
scatter(crd,rrd,5,mygraylines,'.','DisplayName',"River Power Projects")%,'Color',mycolor_riverpower)
bubblelegend("Reservoir volume in MCM")

%% Plot subbasins+newres+ power capacity by symbolsize
figure
h0=imagescnan(catchments_cl);
set(h0, 'AlphaData', .5*(catchments_cl>0))
hold on
tmp1=108+(RD_allSum>0);
tmp1(isnan(RD_allSum))=nan;
h1=imagescnan(tmp1); % all nan vals are set to transparent in this function
colormap([ cmap8; mlabblue;])
cb=colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'TickLabels',basinlabels.basinnames(selbasins),'TickDirection','none','Limits',[101 108] );%'Location','southoutside',
%     % Make the colorbar transparent - Need to do this manually
%     cdata = cb.Face.Texture.CData;
%     cdata(end,:) = uint8(myalpha * cdata(end,:));
%     cb.Face.Texture.ColorType = 'truecoloralpha';
%     cb.Face.Texture.CData = cdata;

bubblechart(crd,rrd,RDPnets_GWh,'k' ,'DisplayName',"GWh indicators",'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.4) %,'MarkerEdgeColor',[1 1 1]*.2)
scatter(crd,rrd,5,mygraylines,'.','DisplayName',"River Power Projects")%,'Color',mycolor_riverpower)
bubblelegend("Energy in GWh per year")

%% Plot Q_wR + power capacity by symbolsize + reservoir size color + opacity=COE but this opacity scaling is not so good
figure
imagescnan(Q_burnR);
hold on
%b=bubblechart(crd,rrd,RDPnets_GWh,ResVols_MCM); %,'MarkerFaceAlpha',COEs) %,'MarkerEdgeColor',[1 1 1]*.2)
s=scatter(crd,rrd,RDPnets_GWh,ResVols_MCM,'filled', 'DisplayName',"River Power Projects"); %,'Color',mycolor_riverpower)
s.AlphaData = 1-rescale(sqrt(COEs), 0.4,1); %log(COEs); % high cost has low alpha/visibility
s.AlphaDataMapping='none'; 
s.MarkerFaceAlpha = 'flat';
%s.MarkerEdgeAlpha = 0.8;
%bubblelegend("Energy in GWh per year")
mycbar("Reservoir volume in MCM")

%% OK:Plot Q_wR + power capacity by symbolsize + COE color + opacity=reservoir size
figure
imagescnan(Q_burnR);
hold on
%b=bubblechart(crd,rrd,RDPnets_GWh,ResVols_MCM); %,'MarkerFaceAlpha',COEs) %,'MarkerEdgeColor',[1 1 1]*.2)
s=scatter(crd,rrd,RDPnets_GWh,COEs,'filled', 'DisplayName',"River Power Projects"); %,'Color',mycolor_riverpower)
s.AlphaData = 1-rescale(ResVols_MCM,0.3,0.95); % high volume has low alpha/visibility
s.AlphaDataMapping='none'; 
s.MarkerFaceAlpha = 'flat';
s.MarkerEdgeColor = 'flat';
set(gca,'ColorScale','log')
mycbar("Cost of production (USD/kWh)")
text(40,1600, "Resevoir volume determines transparency. High volume has high transparency")

%% Plot Q_wR + power capacity by symbolsize + COE color
figure
imagescnan(Q_burnR);
hold on
bubblechart(crd,rrd,RDPnets_GWh,COEs,'MarkerFaceAlpha',0.2) %,'DisplayName',"Reservoir size indicators",'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.4) %,'MarkerEdgeColor',[1 1 1]*.2)
scatter(crd,rrd,5,mygraylines,'.','DisplayName',"River Power Projects")%,'Color',mycolor_riverpower)
bubblelegend("Energy in GWh")
colormap(viridis)
set(gca,'ColorScale','log')
mycbar("Cost of production (USD/kWh)")

%% Plot Q_wR + resvol by symbolsize + COE color
% GOOD: plot catchments overlaid w new reservoirs, symbol size- ResVol
figure
% h0=imagescnan(catchments_cl);
% set(h0, 'AlphaData', .5*(catchments_cl>0))
imagescnan(Q_burnR);
hold on
%scatter(crd,rrd,ResVols_MCM,COEs,'o','filled','MarkerFaceAlpha',0.4,'MarkerEdgeColor',[1 1 1]*.2)
bubblechart(crd,rrd,ResVols_MCM,COEs,'DisplayName',"Reservoir size indicators") %,'MarkerFaceAlpha',0.4,'MarkerEdgeColor',[1 1 1]*.2)

plot(crd,rrd,'.b','DisplayName',"River Power Projects")%,'Color',mycolor_riverpower)
plot(coss(SysID==1),ross(SysID==1),'.','Color', mlaborange,'DisplayName',"Diversion Projects")
mycbar("Cost of production (USD/kWh)")
title("New Reservoirs. Symbol size=Reservior Vol, Color=COE")
subtitle(selrun, 'Interpreter','none')
%colormap(spring)
set(gca,'ColorScale','log')
legend()

%% SCATTER: Plot 5-in-one H vs ResVol symbol size= GWh, color =COE for RP and DP type
figure
bubblechart(ResVols_MCM,OptDHs_m,RDPnets_GWh,COEs,'MarkerFaceAlpha',0.6,'MarkerEdgeColor',[1 1 1]*.2,'DisplayName',"River Power Projects")%,    'MarkerFaceColor',[0 .8 .8]) ;%,'x','Color',copts2(3,:))
hold all
plot(ResVols_MCM,OptDHs_m,'.k')
title('Relation between COE (symbol color), GWh (size), Head (y axis) and ResVol (x axis)')
ylabel('Head in m')
xlabel('Reservoir volume in MCM')
grid on; box on
%legend('Location','Southeast')
mycbar("Cost of production (USD/kWh)")
set(gca,'XScale', 'log') %, 'ColorScale','log','YScale', 'log'
xline(TarbelaVol_MCM,'LineStyle','--','Color','r','LineWidth',1.25,'DisplayName',"Tarbela");
yline(TarbelaHeight_m,'LineStyle','--','Color','r','LineWidth',1.25);
sgtitle(selrun,'Interpreter','none')
legend(["River Power Projects","",  "Values for Tarbela Reservoir"])
bubblelegend("Energy in GWh per year")
