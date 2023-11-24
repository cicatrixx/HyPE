% Prepare cost plots with basin level for 12x2 energy (Tech+sust) X 2 TF
% Prepare cost plots with sub-basin wise distinction of mixed sust case
% cost curve vs demands
clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')
ofmatname=fullfile(rootof,'MainFutScenarios.mat');

load(ofmatname, 'pcsout',...
    'runnames','basindata','catchments_cl')

%% Setup
runnames_cl=strrep(strrep(runnames,'_Tech_Fin','_Technical'),'_Sust_RiskAverse','_Sustainable')';
nscens=length(runnames_cl);
copts=[];
% create my color palette with one shade lighter version for Tech case
for ii=1:12 % 3 scens % create lighter shade for each scenario
    copts=[copts; c_ssp(ii,:);  brighten(c_ssp(ii,:),0);  ];
end
copts=[copts; 0 0 0; 0 0 0]; %gray for historical
%copts=[c_ssp; 0 0 0;c_ssp; 0 0 0;];

%% GOOD: Loop HIST+FUT Cost curve w points for DP and RP separately for - assumes tech and sust scenarios are back to back
myscens_2futs=[3:13*2, 1,2; 27:50,1,2]; %idx list for the two futures w historical at last

for seltf=1:2
    fig=figure('Position', get(0, 'Screensize'),'color','w');
    %subplot(1,2,f);
    c=1;
    mylegend={};
    mylegend_scname={};

    for scn_cnt= myscens_2futs(seltf,:)
        scol=copts(c,:);
        if scn_cnt==1 || scn_cnt==2
            cc(c)=plotCostCurve(pcsout{scn_cnt}.COEAlls, pcsout{scn_cnt}.PnetAlls, pcsout{scn_cnt}.SysIDs,scol,'-',myalpha,"HISTORICAL",0); %
        else
            cc(c)=plotCostCurve(pcsout{scn_cnt}.COEAlls, pcsout{scn_cnt}.PnetAlls, pcsout{scn_cnt}.SysIDs,scol,'-',myalpha,"",0); %
        end
        %    plotCostCurveScatter(pcsout{i}.COEAlls, pcsout{i}.PnetAlls, pcsout{i}.SysIDs,scol,1);
        % Createa alphabetical numbering for legend
        mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
        mylegend_scname= [mylegend_scname(:)' {''} {''} runnames_cl(scn_cnt) ];
        c=c+1; %bcoz i may not go serially
    end
    % Add lines for financial cutoffs
    yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(350,0.02, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)

    xlim([0,450])
    ylim([0.01,100])
    ylabel('Per unit production cost (USD2010 per KWh)',FontWeight='bold')
    xlabel('Basin-wide cumulative energy (TWh per year)',FontWeight='bold')
    set(gca,'YScale','log')
    set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
    %set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
    grid on; box on

    lgd=  legend(cc(1:2:12*2),[repelem([" "],1,8), cornernames],'NumColumns',3,'Orientation','vertical', 'Interpreter','none','Location','best');
    lgd.Title.String=strjoin(sspnames,"   ");
    title(tf_full{seltf})
end

sgtitle("Historical and future cost curves")
% myylim=ylim;

%% FINAL: Sub-basinwise cost curves for all ssps for tech and sust, separated by tfs
for seltf=1:2
    figure
    %s(2)=
    cmap8=cmap8_wong_gray;
    %scn=50; % the Full mixed sust scenario
    addfill=0;
    mylegend={};
    mylegend_basname={};

    for scn_cnt= myscens_2futs(seltf,:)
        c=1;

        for basID= 101:108
            scol=cmap8(c,:);c=c+1; %bcoz i may not go serially
            selbas=pcsout{scn_cnt}.subbasinIDs==basID;
            if rem(scn_cnt,2)
                subplot(1,2,1);hold on % only tech
            else
                subplot(1,2,2);hold on  % only sust
            end

            % differentiate future and hist curves
            if scn_cnt==1 || scn_cnt==2  % for hist add labels and different line type
                plotCostCurveScatter(pcsout{scn_cnt}.COEAlls(selbas), pcsout{scn_cnt}.PnetAlls(selbas), pcsout{scn_cnt}.SysIDs(selbas),"k",addfill,"-.",basindata.basinnames{basID-100});
            else
                plotCostCurveScatter(pcsout{scn_cnt}.COEAlls(selbas), pcsout{scn_cnt}.PnetAlls(selbas), pcsout{scn_cnt}.SysIDs(selbas),scol,addfill,'-');
            end
            mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
            mylegend_basname= [mylegend_basname(:)' {''} {''} [basindata.basinnames(basID-100) ]];
        end
    end

    % Apply formatting
    for sp=1:2
        subplot(1,2,sp)
        % Add lines for financial cutoffs
        yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
        text(18,0.03, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)
        xlabel('Subbasin-wise cumulative energy (TWh per year)',FontWeight='bold')
        ylabel('Per unit production cost (USD2010 per KWh)',FontWeight='bold')
        ylim([0.0100 50*seltf])
        set(gca,'YScale','log')
        set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
        grid on; box on
        title(sprintf("%s potentail ",pottypes2{sp}))
    end
    sgtitle(tf_full{seltf})
end

%% FINAL: Sub-basinwise cost curves for only TECH or SUST, separated by ssps, separated by tfs
myscens_techthensust=[3:2:13*2; % tech tf1
    27:2:50;                    % tech tf2
    4:2:13*2;                   % sust tf1
    28:2:50;                    % sust tf2
    ]; %idx list for tech and sust separated for the two futures
%
myylim=[0.01 100];
scnames=compose("Sub-basin %s cost curve for %s",[pottypes2([1 1 2 2])' tf_full([1 2 1 2])']);
cmap8=cmap8_wong_gray;

for seltf=1:size(myscens_techthensust,1)
    figure
    addfill=0;
    mylegend_basname={};
    sspidx=1;sp=1;
    subplot(1,3,sp);hold on
    for scn_cnt= myscens_techthensust(seltf,:)
        c=1;
        % create new subplot after 4 scens for each ssp
        if sspidx==4
            subplot(1,3,sp);hold on
            sspidx=1;
            sp=sp+1;
        end

        for basID= 101:108
            scol=cmap8(c,:);c=c+1; %bcoz i may not go serially
            selbas=pcsout{scn_cnt}.subbasinIDs==basID;
            plotCostCurveScatter(pcsout{scn_cnt}.COEAlls(selbas), pcsout{scn_cnt}.PnetAlls(selbas), pcsout{scn_cnt}.SysIDs(selbas),scol,addfill,'-');
            mylegend_basname= [mylegend_basname(:)' {''} {''} [basindata.basinnames(basID-100) ]];
        end
        sspidx=sspidx+1;
    end

    % Apply historical and formatting in all subplots
    for sp=1:3
        subplot(1,3,sp)
        for basID= 101:108
            if rem(scn_cnt,2) % odd= tech
                scn_cnt=1;
            else
                scn_cnt=2;
            end

            selbas=pcsout{scn_cnt}.subbasinIDs==basID;
            plotCostCurveScatter(pcsout{scn_cnt}.COEAlls(selbas), pcsout{scn_cnt}.PnetAlls(selbas), pcsout{scn_cnt}.SysIDs(selbas),"k",addfill,"-.",basindata.basinnames{basID-100});
            mylegend_basname= [mylegend_basname(:)' {''} {''} [basindata.basinnames(basID-100) ]];

        end
        % Add lines for financial cutoffs
        yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
        text(18,0.03, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)
        xlabel('Subbasin-wise cumulative energy (TWh per year)',FontWeight='bold')
        ylabel('Per unit production cost (USD2010 per KWh)',FontWeight='bold')
        ylim(myylim)
        set(gca,'YScale','log')
        set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
        grid on; box on
        title(sspnames{sp})
    end
    sgtitle(scnames(seltf))
end

%%  subcatchment map
figure
h=imagescnan(catchments_cl);
colormap(cmap8)
% add names for subbasins
%  text(bas_x(1:8),bas_y(1:8),basindata.basinnames,"HorizontalAlignment","center","FontAngle","italic")


%% Create subbasin wise bubble plot for just plant sizes for mixed Sust
scn_cnt=6; % the Full mixed sust scenario
c=1;
figure
subplot(1,3,1)
hold on
for basID= 101:108
    scol=cmap8(c,:);c=c+1; %bcoz id may not go serially
    selbas=pcsout{scn_cnt}.subbasinIDs==basID;
    yy=sort(pcsout{scn_cnt}.PnetAlls(selbas));
    %plot(sort(pcsout{scn}.PnetAlls(selbas)),'.','DisplayName',basindata.basinnames{basID-100});
    bubblechart(1:length(yy),yy,yy,scol,'DisplayName',basindata.basinnames{basID-100},'MarkerFaceAlpha',0.20);
end
xlabel("CellID")
ylabel("GWh per year")
grid on; box on
legend
blgd = bubblelegend('GWh/yr');

%%
subplot(1,3,2)
boxplot(GWh_yr,subbas,'Symbol','.','Whisker',Inf,'label',basindata.basinnames)%,scol,'DisplayName',basindata.basinnames{basID-100},'MarkerFaceAlpha',0.20);
set(gca, 'YScale', 'log')
xlabel("Subbasin")
ylabel("GWh per year")
grid on; box on

subplot(1,3,3)
boxplot(USD_kWh,subbas,'Symbol','.','Whisker',Inf,'label',basindata.basinnames)%,scol,'DisplayName',basindata.basinnames{basID-100},'MarkerFaceAlpha',0.20);
xlabel("Subbasin")
ylabel("COE in USD2010/kWh")
grid on; box on
yline(0.1)

%% Bin data
% select fine gaps for 0-0.1 and 0.1-1 but can be coarse for 1-20, no vals beyond that. There
% are <3000 project in fin pot. so bin below 0.1 should have a lot of steps
% to sufficiently characterize this.
% When I dont plot individual points, i loose a sense of how big each
% project is at the beginning - lets accumulate both sum of project sizes and num of
% projects
scn_cnt=1;
coebinedges=[0;
    linspace(0.01,0.1,200)';
    logspace(log10(0.11),log10(1),50)';
    logspace(log10(1.1),log10(50),20)'];
for scn= 1:nscens %myscens_techthensust(f,1:4)
    for bi=1:length(coebinedges)-1
        if scn_cnt==1
            name_bin(bi,scn_cnt)=sprintf("%0.5f <COE<= %0.5f",coebinedges(bi:bi+1));
        end

        % BASIN LEVEL
        % find lower<vals<=uppper
        selprj=pcsout{scn_cnt}.COEAlls>coebinedges(bi) & pcsout{scn_cnt}.COEAlls<=coebinedges(bi+1);
        nprjs_bin(bi,scn_cnt)= sum(selprj) ;
        if nprjs_bin(bi,scn_cnt)==0
        pnet_bin(bi,scn_cnt)= nan;
        else
        pnet_bin(bi,scn_cnt)= sum(pcsout{scn_cnt}.PnetAlls(selprj));
        end
        
        %SUBBASIN LEVEL
         for basID= 101:108
            selprj=pcsout{scn_cnt}.COEAlls>coebinedges(bi) & pcsout{scn_cnt}.COEAlls<=coebinedges(bi+1) & pcsout{scn_cnt}.subbasinIDs==basID;
            nprjs_subbasin_bin(bi,scn_cnt)= sum(selprj) ;
        if nprjs_subbasin_bin(bi,scn_cnt)==0
            pnet_subbasin_bin(bi,scn_cnt,basID-100)= nan;
        else
            pnet_subbasin_bin(bi,scn_cnt,basID-100)= sum(pcsout{scn_cnt}.PnetAlls(selprj));
        end
         
         end
    end
scn_cnt=scn_cnt+1;
end
nprjs_cum=cumsum(nprjs_bin);
pnet_TWh_cum=cumsum(pnet_bin,'omitnan')/1000;
pnet_subbasin_TWh_cum=cumsum(pnet_subbasin_bin,'omitnan')/1000;
%compiled=table(name_bin, nprjs_bin, pnet_bin, cumsum(nprjs_bin) cumsum(pnet_bin));

disp("Create cost curve data")
%% try plot one corner
figure;plot(pnet_cum(:,3:2:9),coebinedges(2:end))
colororder(copts)
figure; plotEnvelopeMinMaxMean(pnet_TWh_cum(:,3:2:9), c_ssp_main(selrcp,:),3,coebinedges(2:end)',1)
%% Create indices for rcpcorners, reorder pcsout tech then sust
%for both tf back to back
techidx=3:2:50;
sustidx=4:2:50;

startidx=1:4:25;
endidx=4:4:25;

% Get indices for rcpcorners
for i=1:3*2 % 3sspx2tf
    rcpgrps_tech{i}=techidx(startidx(i):endidx(i)); % contains 3ssps for near fut followed by far fut
    rcpgrps_sust{i}=sustidx(startidx(i):endidx(i)); % contains 3ssps for near fut followed by far fut
end

% techcorners=[3 5 7 9
% 11 13 15 17
% 19 21 23 25];
% 
% sustcorners=[3 5 7 9
% 11 13 15 17
% 19 21 23 25];

% pcsout_hist=pcsout([1 2]);
% % mid fut
% pcsout_techthensust{1}=pcsout([3:2:13*2 4:2:13*2;]);
% % far fut
% pcsout_techthensust{2}=pcsout([27:2:50 28:2:50]);


%% GOOD: Cost curve w bands at basin scale for mid and far future
tfidx=[1:3; 4:6]; %3 ssps; 2 tfs
c_ssp_main=c_ssp(4:4:12,:);
c_ssp_main=[c_ssp_main; c_ssp_main];
figure
for seltf=1:2
    subplot(1,2,seltf)
    for selrcp=tfidx(seltf,:)
        l1(selrcp)= plotEnvelopeMinMaxMean(pnet_TWh_cum(:,rcpgrps_tech{selrcp}), c_ssp_main(selrcp,:),0,coebinedges(2:end)',1);
        plotEnvelopeMinMaxMean(pnet_TWh_cum(:,rcpgrps_sust{selrcp}), c_ssp_main(selrcp,:),0,coebinedges(2:end)',1);
    end
    l3=plot(pnet_TWh_cum(:,1:2),coebinedges(2:end),'k--','LineWidth',2,'DisplayName','Historical');
    % Add lines for financial cutoffs
    yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(350,0.02, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)

    xlim([0,450])
    ylim([0.01,1])
    ylabel('Per unit production cost (USD2010 per KWh)',FontWeight='bold')
    xlabel('Basin-wide cumulative energy (TWh per year)',FontWeight='bold')
    set(gca,'YScale','log')
    set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
    %set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
    grid on; box on

%     lgd.Title.String=strjoin(sspnames,"   ");
    title(tf_full{seltf})
end
    lgd=  legend([l1(1:3) l3(1)],[rcpnames,"Historical"], 'Interpreter','none','Location','best');

%% GOOD: Cost curve w bands at sub-basin scale for one RCP and TF
    figure
 for seltf=1:2 % 2TF
     for selrcp=tfidx(seltf,end) % rcp8.5
        c=1;
        for basID= (101:108)-100
                scol=cmap8_wong_gray(c,:);     
            subplot(1,2,1)
            plotEnvelopeMinMaxMean(pnet_subbasin_TWh_cum(:,rcpgrps_tech{selrcp},basID), scol,0,coebinedges(2:end)',1);
            plot(pnet_subbasin_TWh_cum(:,1,basID),coebinedges(2:end),'--','Color',scol,'LineWidth',1.5,'DisplayName','Historical');

            subplot(1,2,2)
             l1(c)=plotEnvelopeMinMaxMean(pnet_subbasin_TWh_cum(:,rcpgrps_sust{selrcp},basID), scol,0,coebinedges(2:end)',1);
            l3=plot(pnet_subbasin_TWh_cum(:,2,basID),coebinedges(2:end),'--','Color',scol,'LineWidth',1.5,'DisplayName','Historical');
   
            c=c+1; %bcoz i may not go serially
        end
    end
    for pp=1:2
        subplot(1,2,pp)
    % Add lines for financial cutoffs
    yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(350,0.02, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)

    %xlim([0,50])
    ylim([0.01,0.2])
    ylabel('Per unit production cost (USD2010 per KWh)',FontWeight='bold')
    xlabel('Cumulative energy (TWh per year)',FontWeight='bold')
    set(gca,'YScale','log')
    set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
    %set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
    grid on; box on
    title(pottypes2{pp})
    end

    lgd=  legend([l1 l3(1)],[basindata.basinnames{1:8},"Historical"], 'Interpreter','none','Location','best');

%     lgd.Title.String=strjoin(sspnames,"   ");
    sgtitle([tf_full{seltf} " " rcpnames{3}])
end
    
%% FINAL: Basin cost curve + demand curve
% Supply
c_ssp_main=c_ssp(4:4:12,:);
c_ssp_main=[c_ssp_main; c_ssp_main];
figure
subplot(3,1,1:2)  
for seltf=2
    for selrcp=tfidx(seltf,:)
        l1(selrcp)= plotEnvelopeMinMaxMean(pnet_TWh_cum(:,rcpgrps_tech{selrcp}), c_ssp_main(selrcp,:),0,coebinedges(2:end)',1);
        plotEnvelopeMinMaxMean(pnet_TWh_cum(:,rcpgrps_sust{selrcp}), c_ssp_main(selrcp,:),0,coebinedges(2:end)',1);
    end
    l3=plot(pnet_TWh_cum(:,1:2),coebinedges(2:end),'k--','LineWidth',2,'DisplayName','Historical');
    % Add lines for financial cutoffs
    yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(350,0.02, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)

    xlim([0,450])
    ylim([0.01,1])
    ylabel('Per unit production cost (USD2010 per KWh)',FontWeight='bold')
    xlabel('Basin-wide cumulative energy (TWh per year)',FontWeight='bold')
    set(gca,'YScale','log')
    set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
    %set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
    grid on; box on

%     lgd.Title.String=strjoin(sspnames,"   ");
    title(strcat("Potential supply in ", tf_full{seltf}))
end
    lgd=  legend([l1(tfidx(seltf,:)) l3(1)],[rcpnames,"Historical"], 'Interpreter','none','Location','best');

%% Load demand curve
% Load basin pop and energy usage
popWouter=readtable(fullfile(rootf,"data\data_prep\Wouter_PopFuture2023" ,"Hist_FutPop.csv"));
energySecReq_TWh=readtable(fullfile(rootf,"data\data_prep\EnergyDemandSupply" ,"EnergyReqs_Subbasin.xlsx"),'Sheet','SubbasinEnergyinTWh');

% In R files, Swat is in row 1 and Kabul in row 2. for robustness check
% basinname order in excel matches that in matlab
popWouter=popWouter([2 1 3:10],[2 1 3:24]);
energySecReq_TWh=energySecReq_TWh([2 1 3:10],:);

yrsWouter=[2015	2020:10:2080];
woutersspnames=["Prosperous", "BAU", "Downhill"]; %ssp1,2,3
for ii=1:3
% Separate pop by ssp
    popSubbasin(:,:,ii)=popWouter{:, contains(popWouter.Properties.VariableNames, ["baseline",woutersspnames{ii}])};
    popIndus(:,ii)= popSubbasin(9,:,ii)';
    popUIB(:,ii)=sum(popSubbasin(1:8,:,ii))';
% Separate energy by ssp    
    energySecReq4ssp_TWh(:,:,ii) =energySecReq_TWh{1:10, contains(energySecReq_TWh.Properties.VariableNames, ["baseline",woutersspnames{ii}])};
    energyReqUIB(:,ii)=energySecReq4ssp_TWh(9,:,ii)';
    energyReqIndus(:,ii)=energySecReq4ssp_TWh(10,:,ii)';
end

% Load Pak consumption vs projection
energyPak_TWh=readtable(fullfile(rootf,"data\data_prep\EnergyDemandSupply" ,"EnergyReqs_Subbasin.xlsx"),'Sheet','PakDemandForecasts','DataRange','A4:E40','VariableNamesRange','A3:E3');
disp("Loaded pop and energy use files")

%% GOOD: Add demand curve vs supply climate impact vs pf impact
subplot(3,1,3) %figure
l(1)=plot(energyPak_TWh{:,2},energyPak_TWh{:,1},"k-",'LineWidth',1.5); % Historical gen    
for ssel=1:3
    hold all
    l(2)=plot(energyReqIndus(:,ssel),yrsWouter,":b",'LineWidth',1.1);
    l(3)=plot(energyReqUIB(:,ssel),yrsWouter,":m",'LineWidth',1.1);
    l(4)=plot(energyPak_TWh{:,ssel+2},energyPak_TWh{:,1},":c",'LineWidth',1.1); % Future gen projections
    ylim([2000 2085])
    xlim([0 450])
end
grid on
xlabel("TWh per year",FontWeight='bold')
ylabel("Years",FontWeight='bold')
legend(["Historical Electricity Generation: Pakistan", "Projected Electricity Security Req: All Indus", "Projected Electricity Security Req: UIB Indus",  "Projected Electricity Demand: Pakistan"]);
set(gca,'YDir','reverse','XAxisLocation','top')

text(450/2,2090,'Potential Demand in the Future','FontWeight','bold','LineStyle','none','HorizontalAlignment','center');
