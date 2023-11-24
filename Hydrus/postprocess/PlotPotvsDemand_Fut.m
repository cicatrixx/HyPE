% Plots basin level demand-supply analysis time-series
clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')
%c_ssp26=[0 0 0; c_ssp26; c_ssp26]; % Add black for historical

ofmatHPgenTS_basintot=fullfile(rootof,'FutScenarios_HPgen_basintot03'); %w/o extension as it is added later

%% Load HP performance TS
selnewpots=[1 4];
for pot=1:2 % 1=tech, 2=sust
    HP_TS{pot}=load(strcat(ofmatHPgenTS_basintot,"_",newpots4{selnewpots(pot)},".mat"),'hist_annualTS_basin','fut_annualTS_basin',...
        "hist_monTS_basin","fut_annualTS_basin", "fut_monTS_basin", 'LTannual_basin',"LTmon_basin" ); % most of these have data in the form of (pf, timestep, qf)
    % Convert to twh
    HP_TS{pot}.hist_annualTS_basin=HP_TS{pot}.hist_annualTS_basin/1000;
    HP_TS{pot}.fut_annualTS_basin=HP_TS{pot}.fut_annualTS_basin/1000;
    fprintf("Loaded HP annual TS for %s pot \n", newpots4{selnewpots(pot)})
end

load(strcat(ofmatHPgenTS_basintot,"_",newpots4{selnewpots(pot)},".mat"),'cl_pfnames','cl_Qnames')
npfs=length(cl_pfnames);
nqfs=length(cl_Qnames);
futmon_ts1=futmon_ts1';
futmon_ts2=futmon_ts2';

% Add corner name to Q
%cl_Qnames(2:end)=compose("%s: %s",[rcpcornernames rcpcornernames]' ,cl_Qnames(2:end));
cl_Qnames(2:end)=[rcpcornernames rcpcornernames]';

%% Load basin pop and energy usage
load(fullfile(rootf,'data','UI','data','futDemand_pop_energy.mat'),'yrsWouter','futpopSubbasin', 'futpopUIB','futpopIndus','energySecReq_subbasin_TWh','energyReqUIB_TWh' , 'energyReqIndus_TWh', 'energyConsumPak_TWh')
%table(yrsWouter',futpopUIB/1e6,futpopIndus/1e6)
%% GOOD: Matrix std for tech and sust pot
fns = ["std"]; %fieldnames(LTannual_basin);
fnnames=["Standard Deviation"];

tfidx=[1:13; 1 14:25];
figure
myr=1;
minmaxstdmat= table(histrcpcornernames(1:13)','VariableNames',{'pfnames'});

for pot=1:2
    currentpotname=newpots4{selnewpots(pot)};
    for fi=1:length(fns)
        for futtf=1:2
            datatmp=squeeze(HP_TS{pot}.LTannual_basin.(fns(fi)))/1000;
            datain=datatmp(tfidx(futtf,:),tfidx(futtf,:));
            subplot(2,2,myr)
            imagesc(datain)
            xline([1.5, 5.5,9.5])%,'LineWidth',1.5)
            yline([1.5, 5.5,9.5])%,'LineWidth',1.5)

            set(gca,'xtick',1:npfs,'xticklabel',histrcpcornernames,'ytick',1:npfs,...
                'yticklabel',histrcpcornernames,'ticklabelInterpreter','none')%,'xticklabelrotation',90)
            set(gca,'YDir','normal')
            title(strcat(currentpotname, '-',fnnames(fi)))

            colormap(flipud(parula)) %- liked parula coz its brighter
            xlabel(["Discharge datasets"])%; tf_full{futtf}])
            if fi==1
                ylabel([tf_full{futtf};"Portfolios"])
            end
            axis square
            colorbar
            % caxis([4.7 15.6])

            % Find r,c for best pf for each Q
            [~,R_forPmax]=max(datain,[],1);  % returns max for each column in input matrix
            [~,C_forQmax]=max(datain,[],2);
            [~,R_forPmin]=min(datain,[],1);  % returns max for each column in input matrix
            [~,C_forQmin]=min(datain,[],2);

            % save ridx for best and worst pf for full matrix size
            bestpf_idx(futtf,:,fi)=tfidx(futtf,R_forPmax);
            worstpf_idx(futtf,:,fi)=tfidx(futtf,R_forPmin);

            hold all
            l(1)=scatter(1:13,R_forPmax,"ro",'DisplayName',"Highest SD portfolio");
            l(2)=scatter(1:13,R_forPmin,20,"ko",'DisplayName',"Lowest SD portfolio");

            fprintf("Mean of change in fut %0.1f  \n",mean(datain(2:end,:),'all'))
            minmaxstdmat.(strcat(tf_full{futtf},"_", currentpotname))= [min(datain,[],2), max(datain,[],2)];
            myr=myr+1;
        end
    end
end
sgtitle(sprintf("STD for hist+futpf under hist+futQs for potential") );
legend(l,'Location','northoutside')%,'NumColumns',2)
mycbar("Annual energy generation in TWh/year");
oxlsfile=fullfile(rootoffut,"RangeinMatrix.xlsx");

writetable(minmaxstdmat,oxlsfile,'Sheet',strcat("std_Twh"))

%% Monthly supply only climate impact vs pf impact
c_4corner=hsv(4); %vga(6); %
%c_4corner=c_4corner([1,2,6,5],:);
figure
sgtitle(sprintf("MONTHLY production Technical and sustainable rcp corner pfs under the four corner Qs"))

histidx=1;
for pot=1:2
    for ssel=1:3
        subplot(3,1,ssel)
        %subplot(2,3,ssel)

        hold all
        %hist pf under histQ
        li=1;
        l1(li)= plot(histmon_ts,HP_TS{pot}.hist_monTS_basin(histidx,:,1)'/1000,'k-','DisplayName',"Histpf",'LineWidth',1.5);

        qselnf=sspgrps{ssel};
        qselff=sspgrps{ssel+3};

        % plot 4 corner pf under rcpQs for NF and FF
        for cr=1:4
            li=li+1;
            l1(li)= plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_monTS_basin(qselnf(cr),:,qselnf))/1000,c_4corner(cr,:),0,futmon_ts1); %, 'DisplayName',cl_Qnames{qsel})
            plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_monTS_basin(qselff(cr),:,qselff))/1000,c_4corner(cr,:),0,futmon_ts2); %, 'DisplayName',cl_Qnames{qsel})
        end
        % plot hist pf under rcpQs NF and FF
        plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_monTS_basin(histidx,:,qselnf))/1000,[0 0 0],0,futmon_ts1); %, 'DisplayName',cl_Qnames{qsel})
        plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_monTS_basin(histidx,:,qselff))/1000,[0 0 0],0,futmon_ts2); %, 'DisplayName',cl_Qnames{qsel})

        applymyplotformat(sspnames{ssel})
        %ylim([10 490])
        %        xlim([2000, 2095])
        xline([2036, 2066],'LineStyle',':')
        ylabel("TWh per month")
    end
end
lgd=legend(l1,["Historical", cornernames]); %, "All Indus Energy Security Req", "UIB Indus Energy Security Req", "Pakistan consumption", "Pakistan Projected Generation Requirement"],'NumColumns',2);
lgd.Title.String="Porfolios";
%% GOOD: supply only climate impact vs pf impact
c_4corner=hsv(4); %vga(6); %
%c_4corner=c_4corner([1,2,6,5],:);
figure
histidx=1;
sgtitle(sprintf("Technical and sustainable rcp corner pfs under the four corner Qs" ))
for pot=1:2
    for ssel=1:3
        %        subplot(3,1,ssel)
        subplot(2,3,ssel)
        hold all
        %hist pf under histQ
        li=1;
        l1(li)  = plot(histyrs_ts,HP_TS{pot}.hist_annualTS_basin(histidx,:,1)','k-','DisplayName',"Histpf",'LineWidth',1.5);

        qselnf=sspgrps{ssel};
        qselff=sspgrps{ssel+3};

        % plot 4 corner pf, for each corner plot TS under rcpQs for NF and FF
        for cr=1:4
            li=li+1;
            l1(li)= plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_annualTS_basin(qselnf(cr),:,qselnf)),c_4corner(cr,:),0,fyrs_ts1); %, 'DisplayName',cl_Qnames{qsel})
            plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_annualTS_basin(qselff(cr),:,qselff)),c_4corner(cr,:),0,fyrs_ts2); %, 'DisplayName',cl_Qnames{qsel})
        end
        % plot hist pf under rcpQs NF and FF
        plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_annualTS_basin(histidx,:,qselnf)),[0 0 0],0,fyrs_ts1); %, 'DisplayName',cl_Qnames{qsel})
        plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_annualTS_basin(histidx,:,qselff)),[0 0 0],0,fyrs_ts2); %, 'DisplayName',cl_Qnames{qsel})

        applymyplotformat(rcpnames{ssel})
        ylim([10 490])
        xlim([2000, 2095])
        xline([2036, 2066],'LineStyle',':')
        ylabel("TWh per year")
    end
end
lgd=legend(l1,["Historical", cornernames]); %, "All Indus Energy Security Req", "UIB Indus Energy Security Req", "Pakistan consumption", "Pakistan Projected Generation Requirement"],'NumColumns',2);
lgd.Title.String="Porfolios";

%histminmax(:,:,(pf-1)*3+ssel)]

%% Get end val and growthrate stats for ENERGY
for pot=1:2;
    histidx=1;
    for qq=1:nqfs
        % handle hist and fut discharge separately as differnt # of timesteps
        if qq==1
            seldata=HP_TS{pot}.hist_annualTS_basin(:,:,qq);
            selgrowthTS=getPrctChangeTS(seldata')';
            HP_TS{pot}.hist_growthrates_TS(:,:,qq)=selgrowthTS;
        else
            seldata=HP_TS{pot}.fut_annualTS_basin(:,:,qq);
            selgrowthTS=getPrctChangeTS(seldata')';
            HP_TS{pot}.fut_growthrates_TS(:,:,qq)=selgrowthTS;
        end

        % Annual growth stats are of same size for hist and fut so data stored w pfs x qfs format    % pass data w TS in each cols instead of rows and switch it back to
        % preserve previous indexing style
        HP_TS{pot}.growthrates_mean(:,qq)=mean(selgrowthTS,2);
        HP_TS{pot}.growthrates_std(:,qq)=std(selgrowthTS,0,2);
        % take last vals in the timeseries for each pf under current qf
        HP_TS{pot}.endvals(:,qq)=seldata(:,end);
    end
%     % Get min amd max mean val for each pf
%     rangeingrowthrates{pot}=table(tfhistrcpcornernames', min(HP_TS{pot}.growthrates_mean,[],2), max(HP_TS{pot}.growthrates_mean,[],2),...
%         min(HP_TS{pot}.growthrates_std,[],2), max(HP_TS{pot}.growthrates_std,[],2),'VariableNames',["pf" "mean-min" "mean-max" "std-min" "std-max"] )
end
disp("Evaluated growth rates")

%% GOOD: Loop heatmap growthrate of change MF and FF separate from historical pf under historical qf
fns = ["growthrates_mean" "growthrates_std" "endvals" ]; %fieldnames(LTannual_basin);
fnnames=["Mean" "Standard Deviation" "EndofPeriodVals"];
pot=2;
currentpotname=newpots4{selnewpots(pot)};
tfidx=[1:13; 1 14:25];
figure
for fi=1:length(fns)
    for futtf=1:2
        datatmp=squeeze(HP_TS{pot}.(fns(fi)));
        datain=datatmp(tfidx(futtf,:),tfidx(futtf,:));
        subplot(2,3,(futtf-1)*3+fi)
        imagesc(datain)
        xline([1.5, 5.5,9.5])%,'LineWidth',1.5)
        yline([1.5, 5.5,9.5])%,'LineWidth',1.5)

        set(gca,'xtick',1:npfs,'xticklabel',histrcpcornernames,'ytick',1:npfs,...
            'yticklabel',histrcpcornernames,'ticklabelInterpreter','none')%,'xticklabelrotation',90)
        set(gca,'YDir','normal')
        title(fnnames(fi))

        colormap(flipud(parula)) %- liked parula coz its brighter
        xlabel(["Discharge datasets"])%; tf_full{futtf}])
        if fi==1
            ylabel([tf_full{futtf};"Portfolios"])
        end
        axis square
        colorbar
        %caxis([71 140])

        % Find r,c for best pf for each Q
        [~,R_forPmax]=max(datain,[],1);  % returns max for each column in input matrix
        [~,C_forQmax]=max(datain,[],2);
        [~,R_forPmin]=min(datain,[],1);  % returns max for each column in input matrix
        [~,C_forQmin]=min(datain,[],2);

        % save ridx for best and worst pf for full matrix size
        bestpf_idx(futtf,:,fi)=tfidx(futtf,R_forPmax);
        worstpf_idx(futtf,:,fi)=tfidx(futtf,R_forPmin);

        hold all
        l(1)=scatter(1:13,R_forPmax,"ro",'DisplayName',"Best portfolio");
        l(2)=scatter(1:13,R_forPmin,20,"ko",'DisplayName',"Worst portfolio");
        %l(3)=scatter(C_forPmax,1:13,"r^",'DisplayName',"Best discharge scenario for each portfolio");
        %l(4)=scatter(C_forPmin,1:13,"k^",'DisplayName',"Worst discharge scenario for each portfolio");
    end
end
sgtitle(sprintf("Annual growth rate in hist+futpf under hist+futQs for %s potential",currentpotname) );
legend(l,'Location','northoutside')%,'NumColumns',2)
mycbar("% change in energy per year");

%% GOOD: scatter plot of Mean and std growth rates now and in the future
fns = ["growthrates_mean" "growthrates_std" "endvals" ]; %fieldnames(LTannual_basin);
fnnames=["Mean" "Standard Deviation" "EndofPeriodVals"];

figure
potsym=["o" "x"];
colororder(c_ssp26(1:13,:))

for pot=1:2;
    ppp=1;
    for futtf=1:2
        for fi=1:2
            subplot(2,2,ppp)
            datatmp=squeeze(HP_TS{pot}.(fns(fi)));
            datain=datatmp(tfidx(futtf,:),tfidx(futtf,:))';
            l0(pot,:)=plot(datain,potsym(pot));           
            hold on
            % make hist bolder
            l0(pot,1).LineWidth=1.5;
            l0(pot,1).LineStyle=":";
            set(gca,'XTick',1:13,'XTickLabel',histrcpcornernames(1:13))

            %xline([1.5, 5.5,9.5])%,'LineWidth',1.5)
            xlabel("Discharge scenarios")
            ylabel([tf_full{futtf} "% growth per year"])
            applymyplotformat(fnnames(fi))

            minmax(ppp,:,pot)=[min(datain,[],"all"), max(datain,[],"all")];
            ppp=ppp+1;
        end
    end
end
sgtitle("Annual growth rate of historical and future pf under different Qs")

%l=legend([l(2:end) ;l(1)],[repelem({' ',},24/2-4) cornernames, "Historical"],'NumColumns',4);
lg=legend([l0(1,2:end) l0(1,1) l0(2, 1)],[repelem({' ',},24/2-4) cornernames, "Historical-Tech", "Historical-Sust"],'NumColumns',4);
lg.Title.String=strcat(strjoin(rcpnames,'      '), 'Portfolios');

%% Summary of historical only GRs
gr_temporalstats_hist_all=[];%table(0,0,'VariableNames',["minGR", "maxGR"]);
% histpf under histqf
curdata=HP_TS{pot}.hist_growthrates_TS(histidx,:,1);
gr_temporalstats_hist_all=[gr_temporalstats_hist_all;  mean(curdata) min(curdata) max(curdata)];

% histpf under futqf
curdata=squeeze(HP_TS{pot}.fut_growthrates_TS(histidx,:,2:end));
gr_temporalstats_hist_all=[gr_temporalstats_hist_all;  mean(curdata)' min(curdata)' max(curdata)'];

%% GOOD: TS of growth rate Gradient of supply only climate impact vs pf impact
c_4corner=hsv(4);%
figure
histidx=1;
sgtitle(sprintf("Technical and sustainable rcp corner pfs under the four corner Qs" ))
% plot and eval stats for each of the 48*2 fut TS
for pot=1
    gr_temporalstats_nf=[];
    gr_temporalstats_ff=[];
    for ssel=1:3
        myhistdata=HP_TS{pot}.hist_growthrates_TS;
        myfutdata=HP_TS{pot}.fut_growthrates_TS;
        %        subplot(3,1,ssel)
        subplot(2,3,ssel)
        hold all
        %hist pf under histQ
        li=1;
        curdata=myhistdata(histidx,:,1);
        l1(li)= plot(histyrs_ts(2:end),curdata','k-','DisplayName',"Histpf",'LineWidth',1.5);
        %gr_temporalstats_hist=[gr_temporalstats_hist;  mean(curdata) min(curdata) max(curdata)];

        %select two futures TSs
        qselnf=sspgrps{ssel};
        qselff=sspgrps{ssel+3};
        % plot 4 corner pf under rcpQs for NF and FF
        for cr=1:4
            li=li+1;
            %nf
            curdata=squeeze(myfutdata(qselnf(cr),:,qselnf));
            l1(li)= plotEnvelopeMinMaxMean(curdata,c_4corner(cr,:),0,fyrs_ts1(2:end)); %, 'DisplayName',cl_Qnames{qsel})
            gr_temporalstats_nf=[gr_temporalstats_nf;  mean(curdata)' min(curdata)' max(curdata)'];

            %ff
            curdata=squeeze(myfutdata(qselff(cr),:,qselff));
            plotEnvelopeMinMaxMean(curdata,c_4corner(cr,:),0,fyrs_ts2(2:end)); %, 'DisplayName',cl_Qnames{qsel})
            gr_temporalstats_ff=[gr_temporalstats_ff;  mean(curdata)' min(curdata)' max(curdata)'];
        end
        % plot hist pf under rcpQs NF and FF
        plotEnvelopeMinMaxMean(squeeze(myfutdata(histidx,:,qselnf)),[0 0 0],0,fyrs_ts1(2:end)); %, 'DisplayName',cl_Qnames{qsel}
        plotEnvelopeMinMaxMean(squeeze(myfutdata(histidx,:,qselff)),[0 0 0],0,fyrs_ts2(2:end)); %, 'DisplayName',cl_Qnames{qsel}

        applymyplotformat(rcpnames{ssel})
        ylim([-40 60])
        xlim([2000, 2095])
        xline([2036, 2066],'LineStyle',':')
        ylabel("% growth per year")
    end
    gr_temporalstats{pot}=[gr_temporalstats_nf;gr_temporalstats_ff];
end
lgd=legend(l1,["Historical", cornernames]);
lgd.Title.String="Porfolios";

%% FINAL: Compile mean-min-max growth rate for each rcp group
for pot=1:2
    gr_temporalstats_nf=[];
    gr_temporalstats_ff=[];
    gr_temporalstats_hist_nf=[];
    gr_temporalstats_hist_ff=[];

    myhistdata=HP_TS{pot}.hist_growthrates_TS;
    myfutdata=HP_TS{pot}.fut_growthrates_TS;
    %hist pf under histQ
    curdata=myhistdata(histidx,:,1);
    gr_temporalstats_hist_nf=[gr_temporalstats_hist_nf;  mean(curdata) min(curdata) max(curdata)];

    for ssel=1:3
        %select two futures TSs
        qselnf=sspgrps{ssel};
        qselff=sspgrps{ssel+3};

        % add  hist pf under rcpQs NF and FF tp top of pile
        %nf
        % there are 4 corner Qs for NF - take mean of them
        curdata=squeeze(myfutdata(histidx,:,qselnf));
        gr_temporalstats_hist_nf=[ gr_temporalstats_hist_nf; mean(curdata,'all')' min(curdata,[],'all')' max(curdata,[],'all')'];
        %ff
        curdata=squeeze(myfutdata(histidx,:,qselff));
        gr_temporalstats_hist_ff=[ gr_temporalstats_hist_ff; mean(curdata,'all')' min(curdata,[],'all')' max(curdata,[],'all')' ];

        %  4 corner pf under rcpQs for NF and FF
        for cr=1:4
            %nf
            curdata=squeeze(myfutdata(qselnf(cr),:,qselnf));
            gr_temporalstats_nf=[gr_temporalstats_nf;  mean(curdata,'all')' min(curdata,[],'all')' max(curdata,[],'all')'];

            %ff
            curdata=squeeze(myfutdata(qselff(cr),:,qselff));
            gr_temporalstats_ff=[gr_temporalstats_ff;  mean(curdata,'all')' min(curdata,[],'all')' max(curdata,[],'all')'];
        end
    end
    gr_temporalstats_futbyrcp{pot}=[gr_temporalstats_nf;gr_temporalstats_ff];
    gr_temporalstats_histbyrcp{pot}=[gr_temporalstats_hist_nf;gr_temporalstats_hist_ff];
end
disp("Eval min max growth rate")

%% Write growthrates to excel
oxlsfile=fullfile(rootoffut,'GrowthRates.xlsx');
for pot=1:2
    curpotname=newpots4{selnewpots(pot)};
    tblfut{pot}=array2table(gr_temporalstats_futbyrcp{pot},"VariableNames",strcat(curpotname,"_",["Mean" "Min" "Max"]));
    tblfut{pot}.pfnames=tfhistrcpcornernames(2:end)';
    tblhist{pot}=array2table(gr_temporalstats_histbyrcp{pot},"VariableNames",strcat(curpotname,"_",["Mean" "Min" "Max"]));
    tblhist{pot}.Qnames=tfhistrcpnames';
end

writetable([tblfut{1}(:,[4 1:3]) tblfut{2}(:,1:3)],oxlsfile,'Sheet',"gr_fut")
writetable([tblhist{1}(:,[4 1:3]) tblhist{2}(:,1:3)],oxlsfile,'Sheet',"gr_hist")

disp("Written to excel min max growth rate")

%% FINAL: Compile mean-min-max annual growth rate for each rcp group
oxlsname="RangesinAnnualGrowthRates.xlsx";
myhistdatalumped= {HP_TS{1}.hist_growthrates_TS; HP_TS{2}.hist_growthrates_TS};
myfutdatalumped={HP_TS{1}.fut_growthrates_TS; HP_TS{2}.fut_growthrates_TS};
CompileMinMaxByRcpGroups(myhistdatalumped, myfutdatalumped, oxlsname)

%% FINAL: Compile mean-min-max annual energy generated for each rcp group
oxlsname="RangesinAnnualEnergyGen.xlsx";
myhistdatalumped= {HP_TS{1}.hist_annualTS_basin; HP_TS{2}.hist_annualTS_basin};
myfutdatalumped={HP_TS{1}.fut_annualTS_basin; HP_TS{2}.fut_annualTS_basin};
CompileMinMaxByRcpGroups(myhistdatalumped, myfutdatalumped, oxlsname)

% %% FINAL: Compile mean-min-max %change in annual energy generated for each rcp group
% oxlsname="RangesinAnnualEnergyGen_prctChange.xlsx";
% myhistdatalumped= {HP_TS{1}.hist_annualTS_basin; HP_TS{2}.hist_annualTS_basin};
% myfutdatalumped={HP_TS{1}.fut_annualTS_basin; HP_TS{2}.fut_annualTS_basin};
% CompileMinMaxByRcpGroups(myhistdatalumped, myfutdatalumped, oxlsname)

%% FINAL: TS of demand vs supply climate impact vs pf impact
c_4corner=hsv(4); %vga(6); %
%c_4corner=c_4corner([1,2,6,5],:);
figure
histidx=1;
sgtitle(sprintf("Technical and sustainable rcp corner pfs under the four corner Qs"))
for pot=1:2
    for ssel=1:3
        %        subplot(3,1,ssel)
        subplot(2,3,ssel)
        hold all
        %hist pf under histQ
        li=1;
        l1(li)= plot(histyrs_ts,HP_TS{pot}.hist_annualTS_basin(histidx,:,1)','k-','DisplayName',"Histpf",'LineWidth',1.5);

        qselnf=sspgrps{ssel};
        qselff=sspgrps{ssel+3};

        % plot 4 corner pf under rcpQs for NF and FF
        for cr=1:4
            li=li+1;
            l1(li)= plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_annualTS_basin(qselnf(cr),:,qselnf)),c_4corner(cr,:),0,fyrs_ts1); %, 'DisplayName',cl_Qnames{qsel})
            plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_annualTS_basin(qselff(cr),:,qselff)),c_4corner(cr,:),0,fyrs_ts2); %, 'DisplayName',cl_Qnames{qsel})
        end
        % plot hist pf under rcpQs NF and FF
        plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_annualTS_basin(histidx,:,qselnf)),[0 0 0],0,fyrs_ts1); %, 'DisplayName',cl_Qnames{qsel})
        plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_annualTS_basin(histidx,:,qselff)),[0 0 0],0,fyrs_ts2); %, 'DisplayName',cl_Qnames{qsel})


        % add energy demand metrics
        if pot==2
            l1(li+1)=plot(yrsWouter,energyReqIndus_TWh(:,ssel),"--b",'LineWidth',1.5);
            l1(li+2)=plot(yrsWouter,energyReqUIB_TWh(:,ssel),":b",'LineWidth',1.5);
            l1(li+3)=plot(energyConsumPak_TWh{:,1},energyConsumPak_TWh{:,2},"Color",mybrown,'LineWidth',1); % Historical gen
            l1(li+4)=plot(energyConsumPak_TWh{:,1},energyConsumPak_TWh{:,ssel+2},":","Color",mybrown,'LineWidth',1.2); % Future gen projections
        end

        applymyplotformat(sspnames{ssel})
        ylim([10 490])
        xlim([2000, 2095])
        xline([2036, 2066],'LineStyle',':')
        ylabel("TWh per year")
    end
end
lgd=legend(l1,["Historical", cornernames, "All Indus Energy Security Req", "UIB Indus Energy Security Req", "Pakistan consumption", "Pakistan Projected Generation Requirement"],'NumColumns',2);
lgd.Title.String="Porfolios                          Energy needs";

%% FINAL: TS of num ppl that can be supported by annual TS for the pfs in their own qfs
% hist_annual_ts is 25pfs x 40 yrs, fut is 25pfs x 30 yrs x 25 qs
figure
sgtitle("Technical and sustainable rcp corner pfs under the four corner Qs")
for pot=1:2
    for ssel=1:3
        %        subplot(3,1,ssel)
        subplot(2,3,ssel)
        hold all
        %hist pf under histQ
        li=1;
        l12(li)= plot(histyrs_ts,HP_TS{pot}.hist_annualTS_basin(histidx,:,1)'/energyreq_MWhpercapita,'k-','DisplayName',"Histpf",'LineWidth',1.5);

        qselnf=sspgrps{ssel};
        qselff=sspgrps{ssel+3};

        % plot each of the 4 corner pf under the 4cornerQs for NF and FF
        for cr=1:4
            li=li+1;
            %NF
            l12(li)= plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_annualTS_basin(qselnf(cr),:,qselnf))/energyreq_MWhpercapita,c_4corner(cr,:),0,fyrs_ts1); %, 'DisplayName',cl_Qnames{qsel})
            %FF
            plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_annualTS_basin(qselff(cr),:,qselff))/energyreq_MWhpercapita,c_4corner(cr,:),0,fyrs_ts2); %, 'DisplayName',cl_Qnames{qsel})
        end
        % plot hist pf under rcp1Qs NF and FF
        plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_annualTS_basin(histidx,:,qselnf))/energyreq_MWhpercapita,[0 0 0],0,fyrs_ts1); %, 'DisplayName',cl_Qnames{qsel})
        plotEnvelopeMinMaxMean(squeeze(HP_TS{pot}.fut_annualTS_basin(histidx,:,qselff))/energyreq_MWhpercapita,[0 0 0],0,fyrs_ts2); %, 'DisplayName',cl_Qnames{qsel})

        applymyplotformat(sspnames{ssel})
        %ylim([200 800]) % for tech pot
        ylim([0 800])
        xlim([2000, 2095])
        xline([2036, 2066],'LineStyle',':')
        if ssel==1
            ylabel(sprintf("Millions of people that can \n be supplied with %.2f MWh per year", energyreq_MWhpercapita))
        end
        l12(li+1)=plot(yrsWouter,futpopIndus(:,ssel)/10^6,"--b",'LineWidth',1.5);
        l12(li+2)=plot(yrsWouter,futpopUIB(:,ssel)/10^6,":b",'LineWidth',1.5);
    end
end

lgd=legend(l12,["Historical", cornernames, "All Indus Population","Upper Indus Population"],'NumColumns',5,'Orientation','horizontal');
lgd.Title.String="Porfolios";

%% FINAL: Compile mean-min-max annual energy generated for each rcp group
oxlsname="RangesinPopSupported.xlsx";
myhistdatalumped= {HP_TS{1}.hist_annualTS_basin/energyreq_MWhpercapita; HP_TS{2}.hist_annualTS_basin/energyreq_MWhpercapita};
myfutdatalumped={HP_TS{1}.fut_annualTS_basin/energyreq_MWhpercapita; HP_TS{2}.fut_annualTS_basin/energyreq_MWhpercapita};
CompileMinMaxByRcpGroups(myhistdatalumped, myfutdatalumped, oxlsname)