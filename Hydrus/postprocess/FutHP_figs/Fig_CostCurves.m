% Fig_CostCurves
% Created By    : Sanita Dhaubanjar
% Created For	: SustainIndus project WP2 Paper publication
%=========================
% Code for loading pre-processed output files from Hydrus model to create
% final plots used for publication

% Prepare cost plots with basin level for 12x2 energy (Tech+sust) X 2 TF
% Prepare cost plots with sub-basin wise distinction of mixed sust case
% cost curve vs demands

clc
clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

ofmatname=fullfile(rootof,"FutFigs_cl0","FutFigs_cl1","FigData",'MainFutScenarios.mat');
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

%% FINAL: Cost curve w bands at sub-basin scale for one RCP and TF
subplottight = @(m,n,p) subtightplot (m, n, p, [0.05 0.05],[.11 .1],[.09 .02]); %[vgap hgap], marg_h -[low up],marg_w -[l r]

for seltf=2 % 2TF
    figure
    for selrcp=tfidx(seltf,end) % rcp8.5
        c=1;
        for basID= (101:108)-100
            scol=cmap8_wong_gray(c,:);
            %plot technical
            subplottight(1,2,1)
            plotEnvelopeMinMaxMean(pnet_subbasin_TWh_cum(:,rcpgrps_tech{selrcp},basID), scol,0,coebinedges(2:end)',1);
            plot(pnet_subbasin_TWh_cum(:,1,basID),coebinedges(2:end),'--','Color',scol,'LineWidth',1.5,'DisplayName','Historical');
%plot sustainable
            subplottight(1,2,2)
            l1(c)=plotEnvelopeMinMaxMean(pnet_subbasin_TWh_cum(:,rcpgrps_sust{selrcp},basID), scol,0,coebinedges(2:end)',1);
            l3(c)=plot(pnet_subbasin_TWh_cum(:,2,basID),coebinedges(2:end),'--','Color',scol,'LineWidth',1.5,'DisplayName','Historical');

            c=c+1; %bcoz i may not go serially
        end
    end
    for pp=1:2
        subplottight(1,2,pp)
        % Add lines for financial cutoffs
        yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
        text(350,0.02, sprintf("\\it Finanical potential\n<= %0.2f USD per kWh",costlim),'FontSize',12)

        %xlim([0,50])
        ylim([0.01,0.2])
        ylabel('Production cost (USD2010 per KWh)',FontWeight='bold')
        xlabel('Cumulative energy (TWh yr^{-1})',FontWeight='bold')
        set(gca,'YScale','log')
        set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
        %set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
        grid on; box on
        title(pottypes2{pp})
    end
    text(35,0.02, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)

    lgd=  legend([l3 l1],[repmat([""],1,8), basindata.basinnames{1:8}], 'Interpreter','none','Location','best','NumColumns',2);
    lgd.Title.String=["Historical RCP 8.5"];
    sgtitle(strcat(tf_full{seltf}, " ", rcpnames{3}))
end
    
%% FINAL: Basin cost curve + demand curve
% Supply
c_ssp_main=c_ssp(4:4:12,:);
c_ssp_main=[c_ssp_main; c_ssp_main];
figure
subplot(6,1,1:3)  
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
    ylabel('Production cost (USD2010 per KWh)','FontWeight','bold')
    xlabel('Basin-wide cumulative energy (TWh yr^{-1})','FontWeight','bold')
    set(gca,'YScale','log')
    set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
    %set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
    grid on; box on

%     lgd.Title.String=strjoin(sspnames,"   ");
    title(sprintf("Hydropower availability in %s future", tname{seltf}))
end
    lgd1=  legend([l3(1) l1(tfidx(seltf,:)) ],["Historical" rcpnames], 'Interpreter','none','Location','best');

% GOOD: Add demand curve vs supply climate impact vs pf impact
subplot(6,1,4:5); %
hold all
l(1)=plot(energyPak_TWh{:,2},energyPak_TWh{:,1},"k-",'LineWidth',1.5); % Historical gen    
l(2)= plotEnvelopeMinMaxMean(energyReqIndus(:,1:3),[0 0 1],0,yrsWouter,1,":"); % Wouter projections for Indus
l(3)= plotEnvelopeMinMaxMean(energyReqUIB(:,1:3),'magenta',0,yrsWouter,1,":"); % Wouter projections for UIB
selnotnans=~isnan(energyPak_TWh{:,3});
l(4)= plotEnvelopeMinMaxMean(energyPak_TWh{selnotnans,3:5},'cyan',0,energyPak_TWh{selnotnans,1}',1,":"); % Pak Future gen projections
set(gca,'YDir','reverse','XAxisLocation','top')

applymyplotformat('')
ylim([2000 2085])
xlim([0 450])
xticklabels([])
%xlabel(energylabel,'FontWeight','bold')
ylabel("Years",'FontWeight','bold')
lgd2=legend(l,["Historical Electricity Generation: Pakistan",...
 "Projected Electricity Security Req: All of Indus", ...
 "Projected Electricity Security Req: All of upper Indus", ...
 "Projected Electricity Demand: Pakistan"]);

lgd1.Title.String="Energy availability scenarios";
lgd2.Title.String="Energy requirement scenarios";

text(450/2,2090,'Energy requirements in the future','FontWeight','bold','LineStyle','none','HorizontalAlignment','center');
