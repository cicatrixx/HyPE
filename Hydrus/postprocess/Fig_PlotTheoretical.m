% Plot stacked bars for HP size wise theoretical potential with river spacing
% Plots visualized + theoretical in spatial and bar plots for theoretical

clear all
close all
res="500m"; cellsz_m=500;

minQ=0.1; %m3/s
savefigs=0;
plotfigs=0;
%subplot = @(m,n,p) subtightplot (m, n, p, [0.003 0.02]);
addpath('G:\SurfDrive\HPmodel\Hydrus\postprocess\','G:\SurfDrive\GitConnect\Hydrus\functions',genpath('G:\SurfDrive\GitConnect\Hydrus\devFiles\'))
run('myVarNames.m')

%% Load dataset
root="G:/SurfDrive/GitConnect/output/fig_theoretical_fixedDEM_QbasedChannel";
load(sprintf("%s/%s/TheoreticalPot.mat",root,res));
cmap_vispot =linspecer(3);

datapath='G:/SurfDrive/GitConnect/data/UI/data';
data='UI500m_ArcGIS.mat';
fname=fullfile(datapath,data);
load(fname,'catchments','basinlabels', 'channel','outside','dem')
catchments_cl=maskBasin(catchments,~outside);

%% Load literature values
myval_allbasincells=sum(cell_basinEnergy,'all','omitnan'); %GWh
litvals=readtable("G:\PaperData\LitVals_HPpotential.xlsx",'Sheet', 'Clean','Range','A:B');
litvals_UI=[myval_allbasincells;	litvals.GWh(:)]/1000; %in GWh converted to TWh
litsrc= [{'Theoretical in current study'}  litvals.MatlabLabel(:)'];

%% Get sorted version of channel potentials and set <=0 to nan
dorangekm=dorange*cellsz_m/1000; % specified in terms of km
dorangename=strcat(string(dorangekm'),' km');

ncell=numel(cell_basinEnergy);
for i=1:length(dorangekm)
    sorted_Pot(:,i)=sort(reshape(channel_basinEnergy(:,:,i),[ncell,1]));
    sorted_Pot(sorted_Pot(:,i)<=0,i)=nan;
end

%% Get reshaped version
[nr, nc, nz]=size(channel_basinEnergy);
channel_re=reshape(channel_basinEnergy,[nr*nc,nz]);
channel_re(channel_re<=0)=nan;

%% Tabulate total pot and number of sites in each potential class - all 6 classes
tblsum=table(flip(HPclass)','VariableNames',{'HPtype'});
tblcnt=tblsum;
for i = 1:length(dorangename)
    HPtype = discretize(channel_re(:,i),[0 flip(HPsz_greaterthan) Inf],'categorical',flip(HPclass));
    [potsum, cnt] = grpstats(channel_re(:,i),HPtype,{'sum','numel'});
    tblsum = addvars(tblsum,potsum,'NewVariableNames',dorangename(i));
    tblcnt = addvars(tblcnt,cnt,'NewVariableNames',dorangename(i));
    %tbltmp=grpstats(table(channel_re(:,i), HPtype), 'HPtype',{'sum'},'VarNames',{});
end

disp('% of small project (<50MW) in total')
prct_totpot=sum(tblsum{1:4,2:end})./sum(tblsum{:,2:end})*100
prct_count=sum(tblcnt{1:4,2:end})./sum(tblcnt{:,2:end})*100

%% FINAL: Barchart + line plot on distribution of potential - all classes
figure;
colormap(cbrewer2('BrBG',length(tblsum.HPtype)+1))
%subplot(1,2,1)
b=bar(tblsum{:,2:end}'/1e3,'stacked','FaceColor','flat','FaceAlpha',baralpha,'EdgeAlpha',0);
for k = 1:size(tblcnt{:,2:end}',2)
    b(k).CData = k;
end
%legend(tblsum.HPtype)
xticks(1:numel(dorangekm))
xticklabels(dorangename)
grid on
title("Total potential in each size category")
ylabel('Energy in TWh per year)','FontWeight','bold')
hold all

% GOOD: Proportion of small (<50MW) projects
yyaxis right
b1=plot(prct_totpot,".:",'MarkerSize',20, 'DisplayName','Proportion of small projects','LineWidth',2);
grid on
ylabel('Proportion of plants <50 MW as % of total','FontWeight','bold')
xlabel('River segment length','FontWeight','bold')
xticks(1:numel(dorangekm))
xticklabels(dorangename)
xtickangle(45)
%ylim([0 35])

% Add Litvals
yyaxis left
%hold all
yline(litvals_UI(2),"--","LineWidth",1.3,'Color','k'); %Visualized
for i=3:5
    yline(litvals_UI(i),":","LineWidth",2,'Color',.3*[1 1 1]); %Visualized + lit vals
end
text(repmat(dorangekm(1),4,1),litvals_UI(2:5)+25,litsrc(2:5),'FontAngle','italic')

legend([tblsum.HPtype{:} "" "" "" "" "Proportion of small projects"])

%% GOOD: Double Barchart for num projects and for box plot
figure;
% Barchart for # of projects
colormap(cbrewer2('BrBG',length(tblsum.HPtype)+1))

subplot(1,2,1)
b=bar(tblcnt{:,2:end}','stacked','FaceColor','flat','FaceAlpha',baralpha,'EdgeAlpha',0);
legend(tblsum.HPtype)
xticks(1:numel(dorangekm))
xticklabels(dorangename)
ylabel('Number of sites')
grid on
title("Number of sites in each size category")
set(gca, 'YScale', 'log')
% savefig(fullfile(sprintf("%s/%s/TheoreticalPotential_SizeWise.fig",root,res)))

for k = 1:size(tblcnt{:,2:end}',2)
    b(k).CData = k;
end

% Box plot
subplot(1,2,2)
boxplot(channel_re,'Symbol','.','Whisker',Inf)%,'Notch','on')
set(gca, 'YScale', 'log')
ylabel("Theoretical potential at each segment (GWh/yr)")
xticks(1:numel(dorangekm))
xticklabels(dorangename)
grid on
addHPcategoryLines(01,'h',12)
% add mean as well
hold on
 plot(nanmean(channel_re), '*')
 hold off
 
%% Load visualized DB and existing reservoirs
load('G:\SurfDrive\GitConnect\data\UI\data\Existing+UC+All_3.mat', 'existing_dams','existing_reservoirs')
statvals=unique(existing_dams.status);
energy_GWh=existing_dams.GWh;
c_dams_shift=existing_dams.c_dams;
r_dams_shift=existing_dams.r_dams;
setreservoirstatus=zeros(size(existing_reservoirs));
for k=1:height(existing_dams)
    if strcmp(existing_dams.status(k),"Existing")
        setreservoirstatus(existing_reservoirs==existing_dams.id_dams(k)) = 30;
    elseif strcmp(existing_dams.status(k),"Under Construction")
        setreservoirstatus(existing_reservoirs==existing_dams.id_dams(k)) = 15;
    else
        setreservoirstatus(existing_reservoirs==existing_dams.id_dams(k)) = 15;
    end
end

%% Setup for spatial maps
sel_basinEnergy=channel_basinEnergy(:,:,1); %sel the do=1 case
idxo=  find(~isnan(sel_basinEnergy));
potval=sel_basinEnergy(idxo);
minscal=min([potval; energy_GWh ]);
maxscal=max([potval; energy_GWh ]);
%keyboard

[ro,co]= ind2sub(size(sel_basinEnergy),idxo  );
[rch,cch]=find(channel);
nvalid=numel(ro);
inbasin=single(~outside);
inbasin(outside)=nan;

cmap8=cmap8_wong;
selbasins=1:8;

%% GOOD: Plot subbasins + reservoirs + vizualized DB + theoretical potential
figure
h0=imagescnan(catchments_cl);
set(h0, 'AlphaData', .5*(catchments_cl>0))
hold on
tmp=(setreservoirstatus==30)*109;
tmp(setreservoirstatus~=30)=nan;
h1=imagescnan(tmp); % all nan vals are set to transparent in this function
colormap([cmap8; mlabblue;])
cb=colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'TickLabels',basinlabels.basinnames(selbasins),'TickDirection','none','Limits',[101 108] );%'Location','southoutside',
%     % Make the colorbar transparent - Need to do this manually
%     cdata = cb.Face.Texture.CData;
%     cdata(end,:) = uint8(myalpha * cdata(end,:));
%     cb.Face.Texture.ColorType = 'truecoloralpha';
%     cb.Face.Texture.CData = cdata;

hold all
minmarker=1;
maxmarker=1500;
cmap_vispot=[56/255 87/255 35/255; flipud(cbrewer2('Greys',5))]; %green and greys
%cmap = (cbrewer('div','RdYlBu',6));  %viridis(4);

mlaborange=[0.8500 0.3250 0.0980];
mlabblue=[0 0.4470 0.7410];
scatter(co,ro,rescale(potval,minmarker,maxmarker,'InputMax',maxscal,'InputMin',minscal),mlaborange,'filled','MarkerFaceAlpha',0.5,'DisplayName','Theoretical Potential') %,'MarkerFaceColor','#D95319')cmap(1,:)
%scatter(c_dams_shift,r_dams_shift,rescale(energy_GWh,1,1500,'InputMax',maxscal,'InputMin',minscal),[0 0.4470 0.7410]	,'x','LineWidth',1.15,'DisplayName','Visualized potential')

mkrs=["o", "+","x","x"];
c=1;
for i =[1 4 2 3] %Re-order status
    if i==2
        selidx=strcmp(existing_dams.status,statvals(2)) | strcmp(existing_dams.status,statvals(3));%c=4;
    else
        selidx=strcmp(existing_dams.status,statvals(i));
    end %
    scatter(c_dams_shift(selidx),r_dams_shift(selidx),rescale(energy_GWh(selidx),minmarker,maxmarker,'InputMax',maxscal,'InputMin',minscal),cmap_vispot(c,:),mkrs(c),'LineWidth',1.2,'DisplayName',sprintf('Visualized: %s',statvals{i}))
    c=c+1;
end
legend()
title("Theoretical potential compared to visualized potential")

%% Scatter: Plot DEM + vizualized DB + theoretical potential
figure
minmarker=3;
maxmarker=1800;
imagescnan(maskBasin(dem,~outside))
mycbar('Elevation (m)')
%colormap(flipud(cbrewer('seq','Blues',5)))
cdem=colormap("gray");
colormap(cdem(150:end,:))

hold all
%cmap_vispot = flipud(cbrewer('div','PRGn',4));  %viridis(4);
%cmap_vispot=cmap_vispot([1 4 3],:)*0.8;
%cmap = (cbrewer('div','RdYlBu',6));  %viridis(4);

mlaborange=[0.8500 0.3250 0.0980];
scatter(co,ro,rescale(potval,minmarker,maxmarker,'InputMax',maxscal,'InputMin',minscal),mlaborange,'filled','MarkerFaceAlpha',0.3,'DisplayName','Theoretical Potential') %,'MarkerFaceColor','#D95319')cmap(1,:)

c=1;
for i =[1 4 2 3] %Re-order status
    if i==2
        selidx=strcmp(existing_dams.status,statvals(2)) | strcmp(existing_dams.status,statvals(3));%c=4;
    else
        selidx=strcmp(existing_dams.status,statvals(i));
    end
    scatter(c_dams_shift(selidx),r_dams_shift(selidx),rescale(energy_GWh(selidx),minmarker,maxmarker,'InputMax',maxscal,'InputMin',minscal),cmap_vispot(c,:),'o','LineWidth',1.2,'DisplayName',sprintf('Visualized: %s',statvals{i}))

    c=c+1;
end
legend('Color',cdem(length(cdem)-10,:))
%% GOOD: Bubblechart: Plot DEM + vizualized DB + theoretical potential
figure
imagescnan(maskBasin(dem,~outside))
mycbar('Elevation (m)')
%colormap(flipud(cbrewer('seq','Blues',5)))
cdem=colormap("gray");
colormap(cdem(150:end,:))

hold all
%cmap_vispot = flipud(cbrewer('div','PRGn',4));  %viridis(4);
%cmap_vispot=cmap_vispot([1 4 3],:)*0.8;
%cmap = (cbrewer('div','RdYlBu',6));  %viridis(4);

bubblechart(co,ro,potval,mlaborange,'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.1,'DisplayName','Theoretical Potential') %,'MarkerFaceColor','#D95319')cmap(1,:)

c=1;
for i =[1 4 2 3] %Re-order status
    if i==2
        selidx=strcmp(existing_dams.status,statvals(2)) | strcmp(existing_dams.status,statvals(3));%c=4;
    else
        selidx=strcmp(existing_dams.status,statvals(i));
    end
    bubblechart(c_dams_shift(selidx),r_dams_shift(selidx),energy_GWh(selidx),cmap_vispot(c,:),'MarkerFaceAlpha',0,'LineWidth',1.2,'DisplayName',sprintf('Visualized: %s',statvals{i}))

    c=c+1;
end
legend('Color',cdem(length(cdem)-10,:))
bubblesize([.5 80])
title("Theoretical potential compared to visualized potential")
%% GOOD: Bubblechart: Plot existing Reservoir + vizualized DB + theoretical potential
figure
    imagescnan(inbasin*5-(tmp>0))
    colormap(flipud(cbrewer('seq','Blues',7)))
hold all
%cmap_vispot = flipud(cbrewer('div','PRGn',4));  %viridis(4);
%cmap_vispot=cmap_vispot([1 4 3],:)*0.8;
%cmap = (cbrewer('div','RdYlBu',6));  %viridis(4);

bubblechart(co,ro,potval,mlaborange,'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.1,'DisplayName','Theoretical Potential') %,'MarkerFaceColor','#D95319')cmap(1,:)

c=1;
for i =[1 4 2 3] %Re-order status
    if i==2
        selidx=strcmp(existing_dams.status,statvals(2)) | strcmp(existing_dams.status,statvals(3));%c=4;
    else
        selidx=strcmp(existing_dams.status,statvals(i));
    end
    bubblechart(c_dams_shift(selidx),r_dams_shift(selidx),energy_GWh(selidx),cmap_vispot(c,:),'MarkerFaceAlpha',0,'LineWidth',1.2,'DisplayName',sprintf('Visualized: %s',statvals{i}))

    c=c+1;
end
legend('Color',cdem(length(cdem)-10,:))
bubblesize([.5 80])
title("Theoretical potential compared to visualized potential")
