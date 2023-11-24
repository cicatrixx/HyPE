%% Create plots comparing Q_LTmonavgs at diff outlets for observed, historical and simulated future
% Created By    : Sanita Dhaubanjar on 02 Jun 2020
% Created For	: SustaIndus project WP2
%=========================

clc
close all
clear all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

root=pwd; %"G:/GitConnect";
runTS_FDCextraction=0; % load(=0) or run(=1)
datapath=fullfile(root,"data","UI","data");
matout=fullfile(rootoffut, "RQatoutlets_m3s_02.mat");
%xlsfile_summary=fullfile(root,"output","FutRuns","FDC_change.xlsx");

%% Load downscaled 500m basinmask, catchments, outlets
data500m='UI500m_ArcGIS.mat';
load(fullfile(datapath,data500m), 'outlets','outlet_orig','catchments','outside','basinlabels');
outlets=maskBasin(outlets,~outside,nan);
catchments=maskBasin(catchments,~outside,nan);
basin=maskBasin(double(catchments>0),~outside,nan);

%find all outlets within basin
outletID=unique(outlets(outlets>0)); %~isnan(outlet)
nsubbas=length(outletID);
fprintf("Found %d outlets for following catchments:\n",nsubbas)
countUniques(catchments)
disp("Loaded basin catchment and outlet info and observed station data")

%% Get scenario/filenames in correct order
suffix=strcat("_LTavgs_",strrep(tframe,'-','_'),".mat");
prefix="MLAB500m_FutQ_m3day_";
fpaths_Q500m=fullfile(datapath,"Q_LTMavgs",["MLAB500m_40yrClimatology_m3day_R4.mat"; ...
    strcat(prefix, modorder.matlabname,suffix(1)); ...
    strcat(prefix,modorder.matlabname,suffix(2))]);
nQs=length(fpaths_Q500m);

%% Load 500m hist and future Q maps and extract station TS at all outlets
if runTS_FDCextraction
    %% Load observed station data
    obsQ = readmatrix(fullfile(datapath, "obsQ.xlsx"),'Sheet','Obs_climatology','Range','B3:D15');
    stationnames=["Tarbela", "Mangala", "Marala"];
    stationids=[103 104 105];
    %% Load 500m hist and future Q maps
    % fpaths_Q500m=[path2fldrfiles(datapath,"MLAB500m_40*"); path2fldrfiles(datapath,"MLAB500m_FutQ_m3day*")];  % first file listed is historical, rest are future
    Q500m_TS_compiled=zeros(13,nsubbas,nQs,'double'); %3D matrix w LT avgs for 12+1 months  x 12 outlets x 24+1 fut

    for f=1:nQs
        % Get fut name from file name for all except historical run
        if f>1
            tmps=strsplit(fpaths_Q500m{f},'MLAB500m_FutQ_m3day_');
            futname{f}=strrep(tmps{2},'.mat','');
        end
        load(fpaths_Q500m{f},'Q500m_m3day');
        Q_m3s_SPHY=cat(3,Q500m_m3day{:})/(24*60*60); % convert cell to matrix

        % Mask basins to catchments
        for m=1:13
            Q_m3s_SPHY(:,:,m)=maskBasin(Q_m3s_SPHY(:,:,m),~outside);
        end

        % Extract station TS from Q500m
        for pt=1:nsubbas
            % Get 500m TS at outlets
            [r0,c0]=find(outlets==outletID(pt));
            Q500m_TS_compiled(1:13,pt,f) = squeeze(Q_m3s_SPHY(r0,c0,:));
        end
        fprintf("Extracted TS from 500m Qmaps for Q #%d out of %d\n",f,nQs)
    end
    futname{1}="Historical-Sim";

    %% Get percentage change in hist vs fut
    selpts = find(ismember(outletID,stationids)); % get col number for outlet locations that match the 3 stations (old outlet IDs: [3 8 10])
    [hist_pbias, hist_nse] = evalStats(obsQ(1:12,:), Q500m_TS_compiled(1:12,selpts,1));

    for f=2:nQs % #1 is hist
        [Pbias.histobs(f-1,:), ~] = evalStats(obsQ(1:12,:), Q500m_TS_compiled(1:12,selpts,f));
        [Pbias.histsim(f-1,:),~] = evalStats(Q500m_TS_compiled(1:12,selpts,1), Q500m_TS_compiled(1:12,selpts,f));
        Pchange_annual_histsim(f-1,:) = (Q500m_TS_compiled(13,selpts,f)-Q500m_TS_compiled(13,selpts,1))./Q500m_TS_compiled(13,selpts,1)*100;
    end

    %% Evaluate FDC
    obsFDC=sort(obsQ(1:12,:,:),'descend');
    myFDC=sort(Q500m_TS_compiled(1:12,:,:),'descend'); % sort only monthly data
    myFDC_changeprct=(myFDC-myFDC(:,:,1))./myFDC(:,:,1)*100; % 12 mon x 8 subbas x 25 scens
        %% Evaluate design discharge params for selected exceedances
        designQnames=["Q30", "Q40", "Q80"];
        designQ_exceedances=[30, 40, 80];
        designQ_selmonths=round(designQ_exceedances/100*12);
        designQ_selchanges=[];
        
        % rearrange shape of myFDC_changeprct so subbasin is in 3rd dimension
        for selsubbas=1:nsubbas
            myFDC_changeprct_subbasin(1:12,:,selsubbas)=squeeze(myFDC_changeprct(:,selsubbas,:)); %12mons x 25 scens
            %
            designQ_selchanges_subbasin(1:3,:,selsubbas)=myFDC_changeprct_subbasin(designQ_selmonths,:,selsubbas);
            % Stack changes for selected months each subbasin at a time
            designQ_selchanges=[designQ_selchanges; myFDC_changeprct_subbasin(designQ_selmonths,:,selsubbas)];
            %
            if exist("xlsfile_summary")
                writetable(table(myFDC_changeprct_subbasin(:,:,selsubbas)),xlsfile_summary,'Sheet',basinlabels.basinnames{basinlabels.basinIDs==outletID(selsubbas)});
            end
        end
        
        %% Extract and write dchange for specific FDC to excel
        scnames=strrep(extractAfter(fpaths_Q500m,'m3day'),".mat","");
        designQ_selchanges=array2table(designQ_selchanges,'VariableNames',scnames);
        designQ_selchanges.Qname= repmat(designQnames',nsubbas,1);
        designQ_selchanges.basinname=repelem(basinlabels.basinnames(1:8),length(designQ_selmonths));
        xx=designQ_selchanges(:,[26,1:25]);
        if exist("xlsfile_summary")
            writetable(xx,xlsfile_summary,'Sheet','Design_Qs_prct')
        end
        disp("Evaluated and written dchange FDC")


    %% Save
    if exist('matout')
        save(matout,'Q500m_TS_compiled','futname','obsFDC','myFDC','myFDC_changeprct','designQ_selchanges_subbasin','-append')
    end
    disp("Compiled and saved Q_TS data")

else
    load(matout,'Q500m_TS_compiled', 'futname','obsFDC','myFDC','myFDC_changeprct')
    disp("Loaded compiled Q_TS data")
end

%% Bar plot of change in ANNUAL AVG Q subbasin-wise = same as theoretical potential
%barplot is easier to read than scatter. so kept bar
Q500m_changeprct=(Q500m_TS_compiled-Q500m_TS_compiled(:,:,1))./Q500m_TS_compiled(:,:,1)*100;

idatatot=squeeze(Q500m_TS_compiled(13,:,:));
idataprct=squeeze(Q500m_changeprct(13,:,:));
figure;
subplot(3,4,1:4)
bar(idatatot(1:8,2:end),'EdgeAlpha',0)
xline(1:8,'color','black','LineWidth',1.05) % tf splits
ylabel('Annual average discharge (m^3/s)','fontweight','bold')
xticklabels(basinlabels.basinnames)
xtickangle(0)
applymyplotformat('Sub-basin wise annual average discharge',c_ssp)
l=legend([repelem({' ',},24/2-4) cornernames],'NumColumns',3);
l.Title.String=strjoin(tfrcpnames,'      ');

for b=1:8
    subplot(3,4,b+4)
    bar(reshape(idataprct(b,2:end),12,[])','EdgeAlpha',0)
    %bar(reshape(subPot_prctchange(b,2:end),4,[])') this works but cant change color
    ylim([-25, 95])
    applymyplotformat(basinlabels.basinnames(b),c_ssp)
    %xline([4.5:4:24],'LineStyle',':') % ssp splits
    xline([1.5],'color','black','LineWidth',1.05) % tf splits
    xticklabels(strcat(tname,": ", tframe))
end
subplot(3,4,5)
ylabel("% Change in annual average discharge",'fontweight','bold')
%text([1 2],[91 91],strcat(tname,": ", tframe),'HorizontalAlignment','center')%,'fontweight','bold')

%% BAD: Plot bar chart w % change in annual avg discharge for all 8 stations grouped by scenario
figure
bar(squeeze(Q500m_changeprct(13,:,2:24)))
ylabel("% change in annual average Q historical vs fut")
grid on
set(gca,'xtick',1:nQs-1,'xticklabels',cornernames,'TickLabelInterpreter','none','XTickLabelRotation',0) % label w cornernames
%set(gca,'xtick',1:nQs-1,'xticklabels',futname(2:end),'TickLabelInterpreter','none','XTickLabelRotation',90) % label w scennames
% label ssps
text(2.5:4:24,-25*ones(1,6),repmat(modelnames(3:5),1,2),'HorizontalAlignment','center','FontWeight','bold','FontAngle','italic')
% label time frames
text(6.5,95,tframe(1),'FontWeight','bold','FontSize',16,'HorizontalAlignment','center')
text(12.5+6.5,95,tframe(2),'FontWeight','bold','FontSize',16,'HorizontalAlignment','center')
%xline(12.5,'k','LineWidth',2) % for tf
%xline([4.5:4:24],'Color',mygraylines,'LineWidth',1.15)  % for ssp
legend(basinlabels.basinnames(1:8),'Orientation','horizontal',Location='bestoutside')
colororder(cmap8_wong)
%BETTER TO PLOT THIS BY SUBBASIN USING THE THEORYPOT THING



%% GOOD: Plot time series sim hist vs fut - COLORFUL LINES
figure
colororder( brighten(c_ssp_main,0.3))
for selsubbas=1:nsubbas
    subplot(2,4,selsubbas)
    l3=plot(squeeze(Q500m_TS_compiled(:,selsubbas,2:nQs))); %,"DisplayName","Sim: Future")
    hold on
    l2=plot(Q500m_TS_compiled(:,selsubbas,1),"-ok",'LineWidth',1.25);%,"DisplayName","Sim: Historical");%,'Color',[0 0 0]+.5)

    % Plot observed if available
    obsidx=find(outletID(selsubbas)==stationids);
    if obsidx
        l1=plot(obsQ(:,obsidx,1),".k",'LineStyle','-');%,"DisplayName","Obs: Historical");
    end
    applymyplotformat(basinlabels.basinnames{basinlabels.basinIDs==outletID(selsubbas)})
    set(gca,'xtick',1:13,... %,'linestyleorder', {':o','-.s','-d'}
        'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','ANNUAL'})
    xline(12.5)%,'LineStyle',':')
    ylabel("Discharge in m^3/s")%500m Routed
end
sgtitle("Simulated monthly LT avg TS for historical and future")

% Apply labels and formatting to graph
legend([l1 l2 l3(1:8:24)'],modelnames,'Orientation','horizontal') %,'Interpreter','none')

%orient('landscape'); print('G:\GitConnect\output\TheoryPot_WrapUp\Qsimvsobs.pdf','-dpdf','-fillpage')

%% GOOD: Plot future change % time series sim hist vs fut- COLORFUL LINES
figure
colororder( brighten(c_ssp_main,0.3))
for selsubbas=1:nsubbas
    subplot(2,4,selsubbas)
    l3=plot((squeeze(Q500m_TS_compiled(:,selsubbas,2:nQs))-Q500m_TS_compiled(:,selsubbas,1))./Q500m_TS_compiled(:,selsubbas,1)*100); %,"DisplayName","Sim: Future")
    hold on
    %l2=plot(Q500m_TS_future(:,p,1),"-ok");%,"DisplayName","Sim: Historical");%,'Color',[0 0 0]+.5)

    % Plot observed if available
    obsidx=find(outletID(selsubbas)==stationids);
    if obsidx
        %l1=plot(obsQ(:,z,1),".k",'LineStyle','-');%,"DisplayName","Obs: Historical");
    end
    title(basinlabels.basinnames{basinlabels.basinIDs==outletID(selsubbas)})

    box on
    set(gca,'xtick',1:13,... %,'linestyleorder', {':o','-.s','-d'}
        'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','ANNUAL'})
    xline(12.5)%,'LineStyle',':')
    ylim([-100 620])
    ylabel("% Change in Discharge")
end
sgtitle("Change in Simulated monthly LT avg TS")
% Apply labels and formatting to graph
legend([l3(1:8:24)'],modelnames(3:end),'Interpreter','none','Orientation','horizontal')

%orient('landscape'); print('G:\GitConnect\output\TheoryPot_WrapUp\Qsimvsobs.pdf','-dpdf','-fillpage')

%% GOOD: Plot time series sim hist vs fut w BANDS - Mid and Far future in one
figure
sspidx=[2:4:26];  % start idx for each ssp
sspselect=6:-1:1;
% For each point loop through each of the 3x2SSPs and plot 4 corners
for selsubbas=1:nsubbas
    subplot(2,4,selsubbas)
    hold all
    % Plot future simulations plot from reverse order so overlaps are clearer

    for ssp=sspselect
        l(ssp)=plotEnvelopeMinMaxMean(squeeze(Q500m_TS_compiled(:,selsubbas,(sspidx(ssp):sspidx(ssp+1)-1))),c_ssp_main(ssp,:),0);
    end

    % Plot historical simulation
    hold on
    l2=plot(Q500m_TS_compiled(:,selsubbas,1),"-ok");%,"DisplayName","Sim: Historical");%,'Color',[0 0 0]+.5)

    % Plot observed if available
    obsidx=find(outletID(selsubbas)==stationids);
    if obsidx
        l1=plot(obsQ(:,obsidx,1),".k", 'LineStyle', '-');%,"DisplayName","Obs: Historical");
    end
    title(basinlabels.basinnames{basinlabels.basinIDs==outletID(selsubbas)})

    % Apply labels and formatting to graph
    box on
    set(gca,'xtick',1:13,... %,'linestyleorder', {':o','-.s','-d'}
        'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','ANNUAL'})
    xline(12.5,'LineStyle',':')
    ylabel("Discharge in m^3/s")
end
legend([l1 l2 l],["Observed", "Historical", tfrcpnames(sspselect)],'Interpreter','none')
%sgtitle("500m Simulated monthly LT avg TS for historical and future Q")



%% GOOD: Plot FDC sim hist vs fut w BANDS
sspidx=[2:4:26];  % start idx for each ssp
sspselects=[1:3; 4:6];
% For each point loop through each of the 3x2SSPs and plot 4 corners
for f=1:2
    figure(f*100);clf;
    sspselect=sspselects(f,:);

    for selsubbas=1:nsubbas
        subplot(2,4,selsubbas)
        hold all
        ii=1;
        % Plot future simulations plot
        for ssp=sspselect
            l3(ii)=plotEnvelopeMinMaxMean(squeeze(myFDC(:,selsubbas,(sspidx(ssp):sspidx(ssp+1)-1))),c_ssp_main(ssp,:),0);
            scatter(1:12,squeeze(myFDC(:,selsubbas,(sspidx(ssp):sspidx(ssp+1)-1))),5,c_ssp_main(ssp,:),'o');
            ii=ii+1;
        end
        % Plot historical simulation
        hold on
        l2=plot(myFDC(:,selsubbas,1),"-ok");%,"DisplayName","Sim: Historical");%,'Color',[0 0 0]+.5)

        % Plot observed if available
        obsidx=find(outletID(selsubbas)==stationids);
        if obsidx
            l1=plot(obsFDC(:,obsidx,1),".k", 'LineStyle', '-');%,"DisplayName","Obs: Historical");
        end
        applymyplotformat( basinlabels.basinnames{basinlabels.basinIDs==outletID(selsubbas)})

        % Apply labels and formatting to graph
        xline(designQ_selmonths,'LineStyle',':') % main design discharges
        box on
        xlabel('# of months flow exceeded')
        ylabel("Discharge in m^3/s")
        %decided that log scale didnt add more value
        %set(gca,'YScale','log')
        %ylim([10^1 10^4])
    end
    legend([l1 l2 l3],modelnames)

    sgtitle(sprintf("Historical and future FDC for %s",tf_full{f}))
end

%% GOOD: Plot CHANGE IN FDC sim hist - fut w BANDS
sspidx=[2:4:26];  % start idx for each ssp
sspselects=[1:3; 4:6];
% For each point loop through each of the 3x2SSPs and plot 4 corners
for f=1:2
    figure(f*10);clf;
    sspselect=sspselects(f,:);

    for selsubbas=1:nsubbas
        subplot(2,4,selsubbas)
        hold all
        ii=1;
        % Plot future simulations plot
        for ssp=sspselect
            ll(ii)=plotEnvelopeMinMaxMean(squeeze(myFDC_changeprct(:,selsubbas,(sspidx(ssp):sspidx(ssp+1)-1))),c_ssp_main(ssp,:),0);
            scatter(1:12,squeeze(myFDC_changeprct(:,selsubbas,(sspidx(ssp):sspidx(ssp+1)-1))),5,c_ssp_main(ssp,:),'o');
            ii=ii+1;
        end
        applymyplotformat( basinlabels.basinnames{basinlabels.basinIDs==outletID(selsubbas)})
        ylim([-56 206])

        % Apply labels and formatting to graph
        xline(designQ_selmonths,'LineStyle',':') % main design discharges
        box on
        xlabel('# of months flow exceeded')
        ylabel("% change")
    end
    legend(ll,rcpnames,'Interpreter','none')
    sgtitle(sprintf("Change in FDC for %s",tf_full{f}))

    subplot(2,4,1)
    text(designQ_selmonths,[1 1 1]*3000,designQnames,'HorizontalAlignment','center')
end


%% Calculate min/max mean for each scenario for each subbasin
% For total pot
for selsubbas=1:8
    for selq=1:3
        idata=squeeze(designQ_selchanges_subbasin(selq,:,selsubbas))';
        rcpgroups=[0 repelem([1:6],4)]';

        % Save mean in 3d mat with rcp x pottype x subbasin
        idata_byrcp(:,:,selq)=[repelem(idata(1),4)' reshape(idata(2:end),[4  6])];
        idata_mean(:,selq,selsubbas)=groupsummary(idata,rcpgroups,'mean') ;
        idata_neg(:,selq,selsubbas)= groupsummary(idata,rcpgroups,'min') - idata_mean(:,selq,selsubbas);
        idata_pos(:,selq,selsubbas)= groupsummary(idata,rcpgroups,'max')- idata_mean(:,selq,selsubbas);
    end
    subbasdata_byrcp_prct{selsubbas}=idata_byrcp;

end
nrcpscens=size(idata_mean,1);

%% FINAL: Range plot for % changes in sub-basin design discharges- NO theoretical
c_sel= cbrewer2("Set2",3); %cmap8_wong_gray;
xcenter=1:nrcpscens;
xshifted=[xcenter-.15;
    xcenter;
    xcenter+.15    ];
figure
for selsubbas=1:8
    nexttile
    hold all

    for selq=1:3
        ee1(selq)=errorbar(xshifted(selq,:),idata_mean(:,selq,selsubbas), idata_neg(:,selq,selsubbas), idata_pos(:,selq,selsubbas),'_', ...
            'MarkerSize',0.1,'CapSize',0,'LineWidth',2,'Color',c_sel(selq,:),'DisplayName',designQnames{selq});
        %ee2=plot(xshifted(selq,:),idata_mean(:,selq,selsubbas),"_k",'DisplayName','Mean of corners'); %,'Color',mygraylines)
        ee2=plot(xshifted(selq,:),subbasdata_byrcp_prct{selsubbas}(:,:,selq),".",'Color',mygraylines,'DisplayName','4 corners'); %,'Color',mygraylines)

    end
    set(gca,'xtick',1:7,'xticklabel',["Historical" tfrcpnames]) %Xticklength
    xline([1.5 4.5]); %,'LineWidth',1.5)
    applymyplotformat(basinlabels.basinnames(selsubbas ))
    set(gca,'XGrid','off','YMinorGrid','on' )
    xlim([1.5 7.5])
    %add red hline at y=0
    yline(0,'Color','red')
    ylabel("% change")

    if selsubbas<9
        ylim([-50 150])
        % ylim([0.1 1100]) % if theory pot incl
    end

    if selsubbas==9
        set(gca,'YColor',extraaxis ,'XColor',extraaxis );
    end

end
legend([ee1 ee2(1)],'Location','southoutside', 'Orientation', 'horizontal')

sgtitle("   .")


