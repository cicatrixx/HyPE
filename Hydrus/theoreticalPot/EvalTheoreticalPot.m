% Evaluate theoretical potential using cell or channel method for either of
% 15s, 500m or 5km resolution. Channel method is implemented for multiple river
% stretch intervals (do). Considers only Q>=0.1m3/s
% Created By    : Sanita Dhaubanjar
% Created For	: SustaIndus project WP2
%=========================
close all; clc
clear all

% Load inputs based on res
res="500m";rr=2;   cellsz_m=500;
%res="15s"; rr=2;    cellsz_m=450;
%res="15sUI"; rr=2;    cellsz_m=450;
%res="5km";rr=1;    cellsz_m=5000;

if res=="500m"
    datapath='G:/SurfDrive/GitConnect/data/UI/data';
    %data=sprintf('UI%s.mat',res); % for Z, acc, fdir and adir
    % Get annual avg Q
    QMLAB="MLAB500m_40yrClimatology_m3day_R3.mat";
    load(fullfile(datapath,QMLAB),'Q500m_m3day','R500m_m3day');
    Q_m3s=Q500m_m3day{13}/(24*60*60);
    R_m3s=R500m_m3day{13}/(24*60*60);
    
    % Get all else
    data='UI500m_ArcGIS.mat';
    fname=fullfile(datapath,data);
    load(fname,'acc','fdir','dem','flowdist','outside','catchments','basinlabels') %, 'channel') this has older channel based on facc
    load(fullfile(datapath,'channel_Qbased.mat'),'channel') ;

    Z_m=single(dem);
elseif res=="5km"
    datapath='G:/SurfDrive/GitConnect/data/UI/';
    data=sprintf('UI%s.mat',res); % for Z, acc, fdir and adir
    % Get annual avg Q
    QSPHY = 'SPHY5km_40yrClimatology.mat';
    load(fullfile(datapath,QSPHY),'Qm_SPHY_m3sec');
    Q_m3s=Qm_SPHY_m3sec{13};
    % Get all else
    fname=fullfile(datapath,data);
    load(fname,'acc','fdir','dem','flowdist','outside','catchment')
    Z_m=single(dem);
    %outside = ~basinmask;
elseif res=="15s"
    fname='G:/SurfDrive/GitConnect/data/ASIA/Basin/Basin_5.mat';
    load(fname,'acc','fdir','dir','Q','Z','flowdist','outside')
    Z_m=single(Z);
    Q_m3s=double(Q);
    catchment= ~outside;
elseif res=="15sUI"
    fname='G:/SurfDrive/GitConnect/data/ASIA/Basin/Basin_500.mat';
    load(fname,'acc','fdir','dir','Q','Z','flowdist','outside')
    Z_m=single(Z);
    Q_m3s=double(Q);
    catchment= ~outside;
end

% Output file
save2fldr = sprintf("G:/SurfDrive/GitConnect/output/fig_theoretical_fixedDEM_QbasedChannel/%s",res);
if ~isfolder(save2fldr); mkdir(save2fldr); end
outmat='TheoreticalPot.mat';

verbose=1;
finalplots=1;
showplot=0;

% set not basin cells to nan
Z_m(outside)=nan;
acc(outside)=nan;

% minQ threshold
minQ=0.1;%0.1m3/s
% select only Q>=1m3/s
Q_m3sorig=Q_m3s; %select annual avg
Q_m3s(Q_m3s<minQ)=nan;
fprintf("Loaded input datasets for %s\n",res)

%% FINAL: Find channel and plot it + impact of Q>minQ
if res~="500m" %now channel is pre loaded in 500m
    ma = nanmean(acc(:));  % avg of only basin cells
    channel = acc>ma & ~outside; % only cells above mean flow accumulation area and inside basin
end
% Channel cells
[rch,cch]=ind2sub(size(acc), find(channel==1));
% Q>1 and channel cells
[rq,cq]=ind2sub(size(acc), find(channel==1 & Q_m3s>=minQ));

if finalplots %&& res=='500m'
    % Plot flow acc filter
    %     figure;
    %     %subplot(1,2,2)
    %     plot(sort(acc(:)))
    %     set(gca, 'YScale', 'log')
    %     grid on;hold all;
    %     yline(ma,'--r');
    %     %yline(ma/3,':r');
    %     ylabel('Accumlation value')
    %     legend('Acc','Mean acc')%, 'mean acc/3')
    %     title(sprintf('Mean acc selects %0.2f%% of basin cells',sum(channel(:))/sum(~outside(:))*100))
    
    % Plot all channel vs selected channel
    figure
    imagescnan(Z_m);axis image; grid on
    hold all; plot(cch,rch,'o') %plot channel
    hold all; plot(cq,rq,'.') %plot Q>=1
    legend(sprintf('Channel cells: %0.2f%% of basin cells',sum(channel(:))/sum(~outside(:))*100), sprintf('Q>=%0.1f m^3/s$ + channel cells: %0.2f%% of basin cells',minQ, length(rq)/sum(~outside(:))*100))
    title("Selected cells overlaid on Z")
    
    %savefig(fullfile(save2file,'Channel_allVsselect.fig'))
end

%% Absolute Potential: Evaluate runoff weighted elevation
if res=="500m" %now channel is pre loaded
    g           = 9.8;          %Gravitational acceleration (m/s2)
    rho         = 1000;         %Density of water (kg/m3)
    hr_yr      = 8760;          %Hrs in year... (24*365 = 8760 hr/yr)
    eta        = 1;             %Efficiency turbine
    
    minZ=min(Z_m,[],'all');
    
    % Grid based results
    runoff_basinEnergy0 = rho * g * eta .* Z_m.* R_m3s* hr_yr *1e-12; % energy in TWh
    runoff_basinEnergy = rho * g * eta .* (Z_m-minZ).* R_m3s* hr_yr *1e-12; % energy in TWh
    fprintf("Total runoff weighted elevation TWh: %0.2f\n", nansum(runoff_basinEnergy0(~outside)))
    fprintf("Total runoff weighted elevation TWh onlyUIB: %0.2f\n", nansum(runoff_basinEnergy(~outside)))
end

%% Cell by Cell Potential: Get potential - no Q threshold applied
[cell_basinEnergy0, ~]=evalCell2CellTheoreticalPot(Z_m, fdir, Q_m3sorig, 0);
%cellTot_GWh=sum(cell_basinEnergy(cell_basinEnergy(:)>=0)) ;  %in GWh sum of not nan and not negative pot
%sgtitle(sprintf("For %s data",res))
%fprintf("Total discharge weighted head TWh in UIB: %0.2f\n", nansum(cell_basinEnergy0(~outside))/1000)

%% Cell by Cell Potential: Get potential - w Q threshold applied
[cell_basinEnergy, cell_Hgross]=evalCell2CellTheoreticalPot(Z_m, fdir, Q_m3s, 0);
%sgtitle(sprintf("For %s data",res))
fprintf("Total discharge weighted head TWh in UIB: %0.2f\n", nansum(cell_basinEnergy(~outside))/1000)

%% Channel Potential: Get potential at 1 cell spacing
verbose=0;
[ch_basinEnergy, ch_Hgross] = evalChannelTheoreticalPot(Z_m, fdir, Q_m3s, 1, channel, flowdist, verbose);
%sgtitle(sprintf("For %s data",res))
fprintf("Total channel potential in TWh in UIB: %0.2f\n", nansum(ch_basinEnergy(~outside))/1000)

%% Summarize and compare Cell by Cell and Channel Potential for rsi=1 and minQ
fprintf('\n\n\n')
% fprintf("Total runoff weighted elevation based potential in TWh: %0.2f\n", nansum(runoff_basinEnergy(~outside)))
% fprintf("Total discharge weighted head based potential in TWh: %0.2f\n", nansum(cell_basinEnergy0(~outside))/1000)
% fprintf("Total discharge weighted head based potential min Q in TWh: %0.2f\n", nansum(cell_basinEnergy(~outside))/1000)
% fprintf("Total discharge weighted head based potential min Q at d0=1 in TWh: %0.2f\n", nansum(ch_basinEnergy(~outside))/1000)

TWhpot=[nansum(runoff_basinEnergy(~outside))
    nansum(cell_basinEnergy0(~outside))/1000
    nansum(cell_basinEnergy(~outside))/1000
    nansum(ch_basinEnergy(~outside))/1000];
Theorypottype={ 'Runoff weighted elevation only UIB'
                'Discharge weighted head-difference'
                'Discharge weighted head-difference + minQ'
                'Discharge weighted head-difference + minQ + channel only'
    };
TheorySum=table(Theorypottype,TWhpot)
% empty non basin cells in cellpot
cell_basinEnergy_ch=cell_basinEnergy;
cell_basinEnergy_ch(~channel)=nan;
%
if naneq(ch_basinEnergy, cell_basinEnergy_ch)
    disp('Theoretical cell by cell total is captured correctly for channel at do=1')
end

% Gernaat Remaining Technical Potential:	182,908 GWh
chTot_GWh=nansum(cell_basinEnergy(:));
disp("GWh/yr of channel and non-channel pot")
chanEnergy=nansum(cell_basinEnergy(channel)) %/cellTot_GWh*100
nochanEnergy=nansum(cell_basinEnergy(~channel))
chanEnergy_prct=nansum(cell_basinEnergy(channel))/chTot_GWh*100
nochanEnergy_prct=nansum(cell_basinEnergy(~channel))/chTot_GWh*100

%% Channel Potential: Get potential @ 25km spacing
% verbose=1;
% [ch_basinEnergy, ch_Hgross] = evalChannelTheoreticalPot(Z_m, fdir, Q_m3s, 25000/500, channel, flowdist, verbose);
% sgtitle(sprintf("For %s data",res))

%% Loop Channel Potential: Loop do and get potential
%dorange=[1,3,5,10:10:100];
%old DO range
%dorange=ceil([500, 1000, 2000, 3000, 4000, 5000, 10000:10000:100000])./cellsz_m; % oldrange
dorange=ceil([500, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 25000, 40000 :20000:100000]./cellsz_m); % specified in terms of m and converted to cells...have to be multiples of 500/5000

% Preallocate space
channel_basinEnergy=[];
channel_Hgross=[];
channel_GWh=nan(size(dorange));
channel_idxo=[];
channel_idx_nbr=[];
channel_numo=nan(size(dorange));

for di=1:numel(dorange)    %For 15s res, 1 cell is ~500m. In Hydrus:50 (25km);
    % Get channel Potential for current do
    [basinEnergy, Hgross, idxo, idx_nbr]= evalChannelTheoreticalPot(Z_m, fdir, Q_m3s, dorange(di), channel, flowdist, 0);
    
    % Get not nan values from matrix in column vector form
    basinEnergy1 = basinEnergy(~isnan(basinEnergy(:)));
    
    % Archive
    channel_basinEnergy(:,:,di)=basinEnergy;
    channel_Hgross(:,:,di)=Hgross;
    channel_GWh(di)=sum(basinEnergy1(basinEnergy1(:)>=0));  %in GWh sum of not nan and not negative pot
    channel_idxo{di}=idxo;
    channel_idx_nbr{di}=idx_nbr;
    channel_numo(di)=numel(idxo);
end

%% FINAL: Plot channel pot changing w spacing and all other theory pot values
figure
plot(dorange, channel_GWh/1000,".-")
xlabel('do (# of cells)')
ylabel('TWh')
title("Different types of theoretical potential")
grid on
for i=1:length(TWhpot)
    yline(TWhpot(i),"--");%,"LineWidth",1.3,'Color',ccol(2,:));
    text(20,double(TWhpot(i)+20),sprintf('%s: %0.0f TWh/yr', Theorypottype{i},TWhpot(i)))%,'FontWeight','bold')
    
end
% savefig(fullfile(save2fldr,'DiffTypesOfTheoryPot.fig'))

%% Plot channel pot as a % of cell-by-cell
if showplot
    figure
    plot(dorange, channel_GWh/chTot_GWh*100,".-")
    xlabel('do (# of cells)')
    ylabel('% of theoretical energy')
    title(sprintf("Change in channel potential with spacing for %s data",res))
    grid on
end

%% Output: Save estimated potential
if exist('save2fldr')
    fprintf("Saving data to archive for %s\n",res)
    save(fullfile(save2fldr,outmat), 'runoff_basinEnergy','cell_basinEnergy', 'cell_Hgross',...
        'channel_basinEnergy', 'channel_Hgross', 'channel_GWh', ...
        'channel_idxo', 'channel_idx_nbr', 'channel_numo','dorange', 'TheorySum');
    
    %% Output: Save out to geotiff
    savemat2Pantpetiff(fullfile(save2fldr,'TheoreticalPot_Cell.tif'), cell_basinEnergy)
    savemat2Pantpetiff(fullfile(save2fldr,'TheoreticalPot_Channel_1cell.tif'), ch_basinEnergy(:,:,1))
end

%% FINAL: For 500m res Get subbasin-wise potential from rsi=1 channel potential
if finalplots && res=='500m'
    cmap=linspecer(length(basinlabels.basinIDs),'qualitative');
    for st=basinlabels.basinIDs'
        subPot_GWh(st-100) = sum(ch_basinEnergy(ch_basinEnergy > 0 & catchments ==st),'all');
        subPot_GWh_perarea(st-100) = subPot_GWh(st-100)/(sum(catchments ==st,'all')*cellsz_m^2);
    end
    figure
    %total potential
    subplot(2,2,1)
    %bar(subPot_GWh/1000)
    bar(diag(subPot_GWh/1000),'stacked') % for colorful bars
    ylabel('Total energy (TWh/yr)')
    xticklabels(basinlabels.basinnames)
    xtickangle(45)
    applymyplotformat('Sub-basin wise potential at 1 cell spacing',cmap)
    
    % per unit area
    subplot(2,2,2)
    %bar(subPot_GWh/1000)
    bar(diag(subPot_GWh_perarea),'stacked') % for colorful bars
    ylabel('Total energy per unit area (GWh/yr per m^2)')
    xticklabels(basinlabels.basinnames)
    xtickangle(45)
    applymyplotformat('Sub-basin wise potential at 1 cell spacing',cmap)
    
    subplot(2,2,3:4)
    imagescnan(maskBasin(catchments,catchments>0))
    colormap(cmap)
    colorbar('Ticks', basinlabels.basinIDs,'TickLabels',basinlabels.basinnames,'Direction','reverse')
    
end

%% FINAL: Plot range of do, GWh and lit vals
litvals_Indus=[channel_GWh(1) 	296368.32+523812.96	87600 201493	1944876]/1000; %in GWh converted to TWh
litvals_UI=[channel_GWh(1) 	296368.32+523812.96	87600  182908 	 1679935 ]/1000; %in GWh converted to TWh

%litsrc={'Theoretical in current study' 'Visualized for Indus in India (CEA 2020)'	'Visualized for Pakistan (WAPDA 2012)'	'Theoreticalfor Pakistan (Lieftnick et al 1967)'	'Technical for all Indus (Gernaat et al 2017)'	'Theoretical for all Indus (Hoes et al 2017)'};
litsrc={'Theoretical in current study'                                      ;
    'Visualized for Indus in India (CEA 2020) and Pakistan (WAPDA 2012)';
    'Theoretical for Pakistan (Lieftnick et al 1967)'                   ;
    'Technical for Upper Indus at 25km (Gernaat et al 2017)'              ;
    'Theoretical for Upper Indus at 0.22km (Hoes et al 2017)'               };

ccol = linspecer(length(litvals_UI),'qualitative');
do_Gernaat=25; % in km
do_Hoes= 225/1000; % in km (pg 3)
if finalplots
    figure
    plot(dorange*cellsz_m/1000, channel_GWh/1000,"*-b","LineWidth",2)%'Color',ccol(1,:),
    hold all
    %yline(cellTot_GWh,"--");
    yline(litvals_UI(2),"--","LineWidth",1.3,'Color',ccol(2,:)); %Visualized
    yline(litvals_UI(3),"-","LineWidth",2,'Color',ccol(3,:)); %Lieftnick
    yline(litvals_UI(4),":","LineWidth",2,'Color',ccol(4,:)); % Gernaat
    yline(litvals_UI(5),"-","LineWidth",1.3,'Color',ccol(5,:)); %Hoes
    
    %line([0 25], [litvals(5) litvals(5)],'-',"LineWidth",1.3,'Color',ccol(5,:)); % Gernaat
    % Add text labels to ylines and one label for potential at min spacing
    text(repmat(30,4,1),litvals_UI(2:5)*1.1,litsrc(2:5))
    text(dorange(1)*cellsz_m/1000,channel_GWh(1)/1000,sprintf('%0.0f TWh/yr', channel_GWh(1)/1000),'FontWeight','bold','Color','blue') % label max potential IDd
    %Gernaat x=25km
    %Hoes x=225m
    legend(litsrc(1))
    %yline(cellTot_GWh,"--");
    %set(gca, 'YScale', 'log')
    ylabel('Total energy (TWh/yr)')
    xlabel('River spacing (km)')
    title(sprintf("Total theoretical potential at %s and Q>=%0.1fm^3/s",res,minQ))
    grid on
    %yyaxis right
    %plot(dorange*cellsz_m/1000, channel_numo,'x:');
    %ylabel('Number of valid outlets')
    %legend('Channel potential', 'Cell by cell potential','Number of valid outlets')
end
%orient('landscape')
%print('G:/GitConnect/output/TheoryPot_WrapUp/TheoryPot_wSpacing.pdf','-dpdf','-fillpage')

%% Older plots
for olderplots=1
    %% For 500m res Evaluate where Hnegatives lie - using channel potential only for major subbasin
    negH=0;
    if negH && res=='500m'
        cmap = linspecer(length(basinID),'qualitative');
        %% Plot elevation profile and where Hnegatives lie -- for all major subbasin
        selcells=channel& catchments<109;
        idxHneg=ch_Hgross<0;
        figure
        subplot(4,3,1:3)
        scatter(flowdist(selcells)*cellsz_m/1000,ch_Hgross(selcells),10,catchments(selcells),".")
        ylabel('Hloss(m)')
        colormap(cmap)
        grid on
        box on
        
        subplot(4,3,4:12)
        scatter(flowdist(selcells)*cellsz_m/1000,Z_m(selcells),10,catchments(selcells),".")
        hold all
        scatter(flowdist(idxHneg & selcells)*cellsz_m/1000,Z_m(idxHneg & selcells),5,"o",'MarkerEdgeColor',[1 1 1]*.8)%'#A2142F')
        %set(gca, 'YScale', 'log')
        xlabel('Distance from basin outlet (km)')
        ylabel('Z (m)')
        colormap(cmap)
        c=colorbar('Location','SouthOutside','Ticks', basinID,'TickLabels',basinname);%,'Direction','reverse');
        c.Label.String="Sub-basins";
        grid on
        box on
        legend("Channel cells", "Hloss negative cells","Location","Best")
        sgtitle(sprintf("Channel cells in all catchments at %s", res))
        % print('G:/GitConnect/output/TheoryPot_WrapUp/ElevationProfile_allcat.pdf','-dpdf','-fillpage')
        % print('G:/GitConnect/output/TheoryPot_WrapUp/ElevationProfile_allcat_hloss.pdf','-dpdf','-fillpage')
        
        %% Plot elevation profile and where Hnegatives lie -- for swat (102)
        selcells=channel& catchments==102;
        idxHneg=ch_Hgross<0;
        figure
        subplot(4,3,1:3)
        scatter(flowdist(selcells)*cellsz_m/1000,ch_Hgross(selcells),10,catchments(selcells),".")
        ylabel('Hloss(m)')
        colormap(cmap(2,:))
        grid on
        box on
        
        subplot(4,3,4:12)
        scatter(flowdist(selcells)*cellsz_m/1000,Z_m(selcells),10,catchments(selcells),".")
        hold all
        scatter(flowdist(idxHneg & selcells)*cellsz_m/1000,Z_m(idxHneg & selcells),20,"o",'MarkerEdgeColor',[1 1 1]*.8)%'#A2142F')
        xlabel('Distance from basin outlet (km)')
        ylabel('Z (m)')
        colormap(cmap(2,:))
        grid on
        box on
        legend("Channel cells", "Hloss negative cells","Location","Best")
        sgtitle(sprintf("Channel cells in Swat sub-basin at %s", res))
        
        % print('G:/GitConnect/output/TheoryPot_WrapUp/ElevationProfile_swat_hloss.pdf','-dpdf','-fillpage')
        
        %% Plot channel map and where Hnegatives lie
        rr=1;
        win=[1,5]; %# of cells
        
        % Locate Hnegative cells
        Hnegs=cell_Hgross<0;
        [hnr, hnc]=ind2sub(size(cell_Hgross),find(Hnegs));
        newchannel=channel+createBuffer(channel,win(rr))+createBuffer(channel,2*win(rr));
        
        % Count where Hnegs fall
        Hnegs2=newchannel.*Hnegs+Hnegs;
        [C,~,ic] = unique(Hnegs2);
        value_counts = [single(C), accumarray(ic,1)];
        value_prc = value_counts(:,2)/sum(Hnegs(:))*100;
        label_counts={'NaNs', 'In basin', 'In channel+ buffer2', 'In channel+ buffer1', 'In channel', 'Total -ve Hgross'};
        tmp=table(label_counts', [value_counts(:,1); 100],[value_counts(:,2);  sum(Hnegs(:))],[value_prc; 1]);
        mytable = tmp([2,5,4,3,6,1],:);
        mytable.Properties.VariableNames = {'Description','Val','Count','as % of total negatives'}
        prcneg_basin=sum(cell_Hgross(:)<0)/sum(~outside(:))*100
        prcneg_ch=sum(cell_Hgross(channel)<0)/sum(channel(:))*100
        
        figure;imagesc(newchannel+1)
        axis image;grid on
        hold all
        plot(hnc,hnr,'.r','MarkerSize',1)
        title(sprintf("Acc>%0.2f",mean(acc(:))))
        %% Plot Hgross values
        figure
        plot(sort(cell_Hgross(channel)),'.-')
        grid on
        ylabel("Hloss(m) values at channel cells")
        sgtitle(res)
    end
    
    %% Maybe: CDFs
    if showplot
        visdo=[500 5000 10000 50000 100000];
        seldo=find(ismember(dorange, visdo/500));
        figure;hold all
        for i=1:length(seldo)
            tmp=channel_basinEnergy(:,:,seldo(i));
            cdfplot(tmp(tmp>=0))
        end
        
        legend(strcat(num2str(visdo'/1000)," km"))
        set(gca, 'XScale', 'log')
        xlabel("Energy (GWh/yr)")
        %ylabel(" the proportion of the values in x less than or equal to t")
        box on
    end
    
    %% histogram -bad
    if showplot
        visdo=[500 5000 10000 50000 100000];
        seldo=find(ismember(dorange, visdo/500));
        if showplot
            figure;hold all
            for i=1:length(seldo)
                clear hn hx
                tmp=channel_basinEnergy(:,:,seldo(i));
                %histogram(tmp(tmp>0),'DisplayStyle','stairs')
                [hn, hx] = histogram(tmp(tmp>=0)); %// use two-output version of hist to get values
                n_normalized = hn/numel(tmp)/(hx(2)-hx(1)); %// normalize to unit area
                plot(hx, n_normalized);
            end
        end
        legend(strcat(num2str(visdo'/1000)," km"))
        set(gca, 'YScale', 'log', 'XScale', 'log')
    end
    %% Boxplot
    if showplot
        figure
        boxplot(channel_re,'PlotStyle','compact')
    end
    %% MAYBE: histogram
    if showplot
        visdo=[500 5000 10000 50000 100000];
        seldo=find(ismember(dorange, visdo/500));
        [hn, hedges] = histcounts(channel_re);
        hx=filter(1,[0.5 0.5],hedges);      %get midpoint of bins; mean of edges
        
        figure;hold all
        for i=seldo
            [hn, hedges]  = histcounts(channel_re(:,i));
            hx=filter(1,[0.5 0.5],hedges);      %get midpoint of bins; mean of edges
            plot(hx(2:end),hn,'.-')
        end
        %set(gca, 'YScale', 'log')
        %set(gca,'XScale', 'log')
        legend(strcat(num2str(visdo'/1000)," km"))
        applymyplotformat()
        xlabel("Energy (GWh/yr)")
        ylabel("Frequency")
        
    end
end
%%
disp("***************************************EOF***************************************")