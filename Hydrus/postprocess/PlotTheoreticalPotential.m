% Plot theoretical potential changes with river spacing
% Plots visualized + theoretical in spatial and bar plots for theoretical

clear all
close all
res="500m"; cellsz_m=500;

minQ=0.1; %m3/s
savefigs=0;
plotfigs=0;
finalfigs=01;
%subplot = @(m,n,p) subtightplot (m, n, p, [0.003 0.02]);
addpath('G:\SurfDrive\GitConnect\Hydrus\functions',genpath('G:\SurfDrive\GitConnect\Hydrus\devFiles\'))

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

%% Load literature values
myval_allbasincells=sum(cell_basinEnergy,'all','omitnan'); %GWh
litvals=readtable("G:\PaperData\LitVals_HPpotential.xlsx",'Sheet', 'Clean','Range','A:B');
litvals_UI=[myval_allbasincells;	litvals.GWh(:)]/1000; %in GWh converted to TWh
litsrc= [{'Theoretical in current study'}  litvals.MatlabLabel(:)'];

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

%% GOOD FIGS
if finalfigs
    %% FINAL: Plot DEM + vizualized DB + theoretical potential
    figure
    minmarker=3;
    maxmarker=1800;
    imagescnan(maskBasin(dem,~outside))
    mycbar('Elevation (m)')
    %colormap(flipud(cbrewer('seq','Blues',5)))
    cdem=colormap("gray");
    colormap(cdem(150:end,:))

    hold all
    cmap_vispot = flipud(cbrewer('div','PRGn',4));  %viridis(4);
    cmap_vispot=cmap_vispot([1 4 3],:)*0.8;
    %cmap = (cbrewer('div','RdYlBu',6));  %viridis(4);

    mlaborange=[0.8500 0.3250 0.0980];
    scatter(co,ro,rescale(potval,minmarker,maxmarker,'InputMax',maxscal,'InputMin',minscal),mlaborange,'filled','MarkerFaceAlpha',0.3,'DisplayName','Theoretical Potential') %,'MarkerFaceColor','#D95319')cmap(1,:)
    %scatter(c_dams_shift,r_dams_shift,rescale(energy_GWh,1,1500,'InputMax',maxscal,'InputMin',minscal),[0 0.4470 0.7410]	,'x','LineWidth',1.15,'DisplayName','Visualized potential')

    c=1;
    for i =[1 4 2] %Re-order status
        if i==2
            selidx=strcmp(existing_dams.status,statvals(2)) | strcmp(existing_dams.status,statvals(3));%c=4;
        else
            selidx=strcmp(existing_dams.status,statvals(i));
        end
        scatter(c_dams_shift(selidx),r_dams_shift(selidx),rescale(energy_GWh(selidx),minmarker,maxmarker,'InputMax',maxscal,'InputMin',minscal),cmap_vispot(c,:),'o','LineWidth',1.2,'DisplayName',sprintf('Visualized: %s',statvals{i}))
        c=c+1;
    end
    legend('Color',cdem(length(cdem)-10,:))

    %% GOOD: Plot basinreservoirs + vizualized DB + theoretical potential
    figure
    imagescnan(inbasin*5-setreservoirstatus)
    colormap(flipud(cbrewer('seq','Blues',7)))
    %imagescnan(inbasin)
    %colormap(.95*[1 1 1])
    hold all
    scatter(co,ro,rescale(potval,1,1500,'InputMax',maxscal,'InputMin',minscal),[0.8500 0.3250 0.0980],'filled','MarkerFaceAlpha',0.3,'DisplayName','Theoretical Potential') %,'MarkerFaceColor','#D95319')cmap(1,:)
    % scatter(c_dams_shift,r_dams_shift,rescale(energy_GWh,1,1500,'InputMax',maxscal,'InputMin',minscal),[0 0.4470 0.7410]	,'x','LineWidth',1.15,'DisplayName','Visualized potential')

    cmap_vispot = cbrewer('qual','Dark2',8);  %viridis(4);
    cmap_vispot(2,:)=[0 0 0];
    c=1;
    for i =[1 4 2]
        if i==2
            selidx=strcmp(existing_dams.status,statvals(2)) | strcmp(existing_dams.status,statvals(3));
        else
            selidx=strcmp(existing_dams.status,statvals(i));
        end
        scatter(c_dams_shift(selidx),r_dams_shift(selidx),rescale(energy_GWh(selidx),1,1500,'InputMax',maxscal,'InputMin',minscal),cmap_vispot(c,:),'o','LineWidth',1.2,'DisplayName',sprintf('Visualized: %s',statvals{i}))
        c=c+1;
    end
    legend
    title("Theoretical potential compared to visualized potential and existing reservoirs")

    %% FINAL: Barchart + line plot on distribution of potential - all classes
    figure;
    colormap(cbrewer2('BrBG',length(tblsum.HPtype)))
    %subplot(1,2,1)
    b=bar(tblsum{:,2:end}'/1e3,'stacked','FaceColor','flat','FaceAlpha',baralpha);
    for k = 1:size(tblcnt{:,2:end}',2)
        b(k).CData = k;
    end
    %legend(tblsum.HPtype)
    xticks(1:numel(dorangekm))
    xticklabels(dorangename)
    grid on
    title("Total potential in each size category")
    ylabel('Total energy (TWh/yr)')
    hold all

    % GOOD: Proportion of small (<50MW) projects
    yyaxis right
    plot(prct_totpot,".:",'MarkerSize',20, 'DisplayName','Proportion of small projects')
    grid on
    %legend('In terms of total potential','In terms of number of sites')
    ylabel('Number of small project as % of total')
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

    legend([tblsum.HPtype{:} "" "" "" "" "" "Proportion of small projects"])
end
%% TESTFIGS
if plotfigs
    %% GOOD: Plot spatial map of theory pot w subbasins for do=1
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

    %
    figure;
    selbasins=1:8;
    imagescnan(inbasin)
    colorbar('Ticks', basinlabels.basinIDs(selbasins),'TickLabels',basinlabels.basinnames(selbasins),'Direction','reverse')

    hold on
    scatter(co,ro,rescale(potval,1,1500),'filled','MarkerFaceAlpha',0.4, 'DisplayName',sprintf("HP: %0.3g- %0.0f GWh",min(potval),max(potval)))%,'MarkerEdgeColor','blue')
    %plot(cch,rch,'.b','markersize',.05,'DisplayName',"Channel") %plot channel
    legend
    title(sprintf("Channel Potenial: %0.2f TWh at %d sites",sum(potval,'omitnan')/1000,nvalid))
    %h=colorbar('Location','eastoutside');
    %h.Label.String='MW';

    %% TRY: Overlay subbasins and existing reservoirs by setting alpha
    % does not work to overlay DEM and subbasins as the colormap for the
    % two have to be the same so subbasin colors are not seens w the wide
    % elevation data range
    % does not work if h1 matrix has vals other than 0 and 1, the colormap
    % gets thrown off
    cmap8=cmap8_wong;
    figure
    h0=imagescnan(catchments_cl);
    %colormap(cmap8)
    set(h0, 'AlphaData', myalpha*(catchments_cl>0))
    hold on
    h1=imagescnan(setreservoirstatus>0);
    set(h1, 'AlphaData', setreservoirstatus>0)
    colormap([mlabblue;cmap8; ])

    %% TRY: Overlay subbasins and existing reservoirs using imagescnan
    figure
    h0=imagescnan(catchments_cl);
    set(h0, 'AlphaData', myalpha*(catchments_cl>0))
    hold on
    tmp=(setreservoirstatus==30)*109;
    tmp(setreservoirstatus~=30)=nan;
    h1=imagescnan(tmp); % all nan vals are set to transparent in this function
    colormap([cmap8; mlabblue;])
    cb=colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'TickLabels',basinlabels.basinnames(selbasins),'TickDirection','none','Limits',[101 108] );%'Location','southoutside',

    %% TRY: Force existing reservoirs and subbasins to be summed to be displayed together
    cat_w_res=catchments_cl;
    cat_w_res(setreservoirstatus==30)=109;
    figure
    h=imagescnan(cat_w_res);
    colormap([cmap8; mlabblue])
    set(h, 'AlphaData', myalpha*(cat_w_res~=109 & ~outside ) + (cat_w_res==109)) %set alpha of reservoirs as 1 and not outside, not reservoir as 0.6
    cb=colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'TickLabels',basinlabels.basinnames(selbasins),'TickDirection','none','Limits',[101 108] );%'Location','southoutside',

    %% GOOD: Plot subbasins + vizualized DB + theoretical potential
    selbasins=1:8;
    figure
    cmap8=cmap8_wong;
    minmarker=2;
    maxmarker=1800;
    h=imagescnan(catchments_cl);
    colormap(cmap8)
    set(h, 'AlphaData', myalpha*(catchments_cl>0))
    cb=colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'TickLabels',basinlabels.basinnames(selbasins),'TickDirection','none');%'Location','southoutside',
    %     % Make the colorbar transparent - Need to do this manually
    %     cdata = cb.Face.Texture.CData;
    %     cdata(end,:) = uint8(myalpha * cdata(end,:));
    %     cb.Face.Texture.ColorType = 'truecoloralpha';
    %     cb.Face.Texture.CData = cdata;

    hold all
    cmap_vispot = flipud(cbrewer('div','BrBG',4));  %viridis(4);
    cmap_vispot=[ 0 0 0; 0.3 0.6740 0.2;1 0 1];
    cmap_vispot=[0 0 0;0 0 0;0 0 0]; %all black
    cmap_vispot=[56/255 87/255 35/255; flipud(brighten(cbrewer2('Greys',3),-.8))];
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


    %% GOOD: Plot subbasins + reservoirs + vizualized DB + theoretical potential
    selbasins=1:8;
    cat_w_res=catchments_cl;
    cat_w_res(setreservoirstatus==30)=109;
    figure
    h=imagescnan(cat_w_res);
    colormap([cmap8; mlabblue])
    set(h, 'AlphaData', myalpha*(cat_w_res~=109 & ~outside ) + (cat_w_res==109)) %set alpha of reservoirs as 1 and not outside, not reservoir as 0.6
    cb=colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'TickLabels',basinlabels.basinnames(selbasins),'TickDirection','none','Limits',[101 108] );%'Location','southoutside',
    %     % Make the colorbar transparent - Need to do this manually
    %     cdata = cb.Face.Texture.CData;
    %     cdata(end,:) = uint8(myalpha * cdata(end,:));
    %     cb.Face.Texture.ColorType = 'truecoloralpha';
    %     cb.Face.Texture.CData = cdata;

    hold all
    minmarker=2;
    maxmarker=1800;

    cmap_vispot = flipud(cbrewer('div','BrBG',4));  %viridis(4);
    cmap_vispot=[ 0 0 0; 0.3 0.6740 0.2;1 0 1];
    cmap_vispot=[0 0 0;0 0 0;0 0 0]; %all black
    cmap_vispot=[56/255 87/255 35/255; flipud(brighten(cbrewer2('Greys',3),-.8))]; % Green and grays
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

    %% GOOD: Plot scatter points for range of do, TWh and lit vals
    % Load literature values
    myval_allbasincells=sum(cell_basinEnergy,'all','omitnan'); %GWh
    litvals=readtable("G:\PaperData\LitVals_HPpotential.xlsx",'Sheet', 'Clean','Range','A:B');
    litvals_UI=[myval_allbasincells;	litvals.GWh([4 1 3 2])]/1000; %in GWh converted to TWh
    litsrc= [{'Theoretical in current study'}  litvals.MatlabLabel([4 1 3 2])'];

    ccol = linspecer(length(litvals_UI),'qualitative');
    %figure
    subplot(1,2,1)
    plot(dorangekm, channel_GWh/1000,"*-","LineWidth",2,'Color',.3*[1 1 1]);
    hold all
    %yline(litvals_UI(1),"--","LineWidth",1.3,'Color',ccol(1,:)); %This study basin wide sum
    yline(litvals_UI(2),"--","LineWidth",1.3,'Color',ccol(2,:)); %Visualized
    yline(litvals_UI(3),":","LineWidth",1.5,'Color',"k"); %Lieftnick
    yline(litvals_UI(4),":","LineWidth",1.5,'Color',"k"); % Gernaat
    yline(litvals_UI(5),":","LineWidth",1.5,'Color',"k"); %Hoes

    % Add text labels to ylines and one label for potential at min spacing
    text(repmat(dorangekm(5),4,1),litvals_UI(2:5)*1.05,litsrc(2:5),'FontAngle','italic')

    % Add 3 labels for potential at min spacing
    text(dorangekm(1)*1.05,channel_GWh(1)/1000,sprintf('\\leftarrow %0.0f TWh/yr @ 500m', channel_GWh(1)/1000),'FontWeight','bold','Color',.3*[1 1 1]) % label max potential IDd
    % add 4km label
    text(dorangekm(5)*1.05,channel_GWh(5)/1000,sprintf('\\leftarrow %0.0f TWh/yr @ 4km', channel_GWh(5)/1000),'FontWeight','bold','Color',.3*[1 1 1]) % label max potential IDd
    % add 2km label
    text(dorangekm(3)*1.05,channel_GWh(3)/1000,sprintf('\\leftarrow %0.0f TWh/yr @ 2km', channel_GWh(3)/1000),'FontWeight','bold','Color',.3*[1 1 1]) % label max potential IDd

    legend(litsrc{1},"Visualized by basin countries","Other estimates in literature")
    %set(gca, 'YScale', 'log')
    ylabel('Total energy (TWh/yr)')
    xlabel('River segment length (km)')
    title(sprintf("Total theoretical potential at %s and Q>=%0.1fm^3/s",res,minQ))
    grid on

    %% GOOD: Create boxplots for HP size distn varying w do
    %figure;
    subplot(1,2,2)
    boxplot(channel_re,'Symbol','.')%,'Notch','on')

    % add lines and labels for different sizes of HP -- from Siddiqui
    HPclass=["Mega (>1000 MW)", "Large (500-1000 MW)", "Medium (50-500 MW)", "Small (5-50 MW)", "Mini (0.15-5 MW)", "Micro (0.005-0.15 MW)", "Pico (<0.005 MW)"];
    HPsz_greaterthan=[1000	500	50	10	5	0.005]/1000*365*24; %MW converted to GWh assuming year round production

    % from Hoes
    %HPclass=["Large (>50MW)", "Small (5-50 MW)", "Mini (0.15-5 MW)", "Micro (0.005-0.15 MW)", "Pico (<0.005 MW)"];
    %HPsz_greaterthan=[50, 5, 0.15, 0.005]*1000*365*24; %MW converted to kWh
    hold on
    for i=1:length(HPsz_greaterthan)
        yline(HPsz_greaterthan(i),'-.k');
    end

    text(repmat(size(channel_re,2),length(HPclass),1),[HPsz_greaterthan*5 HPsz_greaterthan(end)*.9] ,HPclass)
    xticklabels(strcat(string(dorangekm'),'km'))
    xlabel("River segment length")
    ylabel("Theoretical potential at each segment (GWh)")
    title('Distribution of Project Sizes in Basin Channels')
    set(gca, 'YScale', 'log')
    grid on
    % savefig(fullfile(sprintf("%s/%s/TheoreticalPotential_channel.fig",root,res)))

    %% GOOD: Double Barchart for distribution of potential
    figure;
    colormap(magma)
    subplot(1,2,1)
    b=bar(tblsum{:,2:end}'/1e3,'stacked','FaceColor','flat','FaceAlpha',baralpha);
    for k = 1:size(tblcnt{:,2:end}',2)
        b(k).CData = k;
    end
    %legend(tblsum.HPtype)
    xticks(1:numel(dorangekm))
    xticklabels(dorangename)
    grid on
    title("Total potential in each size category")
    ylabel('Total energy (TWh/yr)')

    % Add Litvals
    hold all
    yline(litvals_UI(1),"--","LineWidth",1.3,'Color','r'); %Basin all cells
    for i=2:5
        yline(litvals_UI(i),":","LineWidth",2,'Color',.3*[1 1 1]); %Visualized + lit vals
    end
    % Add text labels to ylines and one label for potential at min spacing
    text(repmat(dorangekm(5),5,1),litvals_UI(1:5)*1.05,litsrc(1:5),'FontAngle','italic')

    % Add text w arrow
    % text(dorangekm(1)*1.05,channel_GWh(1)/1000,sprintf('\\leftarrow %0.0f TWh/yr @ 500m', channel_GWh(1)/1000),'FontWeight','bold','Color',.3*[1 1 1]) % label max potential IDd
    % % add 4km label
    % text(dorangekm(5)*1.05,channel_GWh(5)/1000,sprintf('\\leftarrow %0.0f TWh/yr @ 4km', channel_GWh(5)/1000),'FontWeight','bold','Color',.3*[1 1 1]) % label max potential IDd
    % % add 2km label
    % text(dorangekm(3)*1.05,channel_GWh(3)/1000,sprintf('\\leftarrow %0.0f TWh/yr @ 2km', channel_GWh(3)/1000),'FontWeight','bold','Color',.3*[1 1 1]) % label max potential IDd

    % Barchart for # of projects
    subplot(1,2,2)
    b=bar(tblcnt{:,2:end}','stacked','FaceColor','flat','FaceAlpha',baralpha);
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

end