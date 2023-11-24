% Fig_AnnualPerformance_Matrix_TS.m
% Created By    : Sanita Dhaubanjar
% Created For	: SustainIndus project WP2 Paper publication
%=========================
% Code for loading pre-processed output files from Hydrus model to create
% final plots used for publication
% Creates annual performance plots - long-term matrices and TS

clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

subplottight = @(m,n,p) subtightplot (m, n, p, [0.08 0.025],[.08 .08],[.08 .01]); %[vgap hgap], marg_h -[low up],marg_w -[l r]

%% Load HP performance TS
matname=fullfile(rootof,'FutScenarios_HPgen_basintot03'); %w/o extension as it is added later
selnewpots=[1 4];
for pot=1:2 % 1=tech, 2=sust
    HP_TS{pot}=load(strcat(matname,"_",newpots4{selnewpots(pot)},".mat"),'hist_annualTS_basin','fut_annualTS_basin',...
        "hist_monTS_basin","fut_annualTS_basin", "fut_monTS_basin", 'LTannual_basin',"LTmon_basin" ); % most of these have data in the form of (pf, timestep, qf)
    % Convert to twh
    HP_TS{pot}.hist_annualTS_basin=HP_TS{pot}.hist_annualTS_basin/1000;
    HP_TS{pot}.fut_annualTS_basin=HP_TS{pot}.fut_annualTS_basin/1000;
    fprintf("Loaded HP annual TS for %s pot \n", newpots4{selnewpots(pot)})
end

futmon_ts1=futmon_ts1';
futmon_ts2=futmon_ts2';

%% Load basin pop and energy usage
load(fullfile(rootf,'data','UI','data','futDemand_pop_energy.mat'),'yrsWouter','futpopSubbasin', 'futpopUIB','futpopIndus','energySecReq_subbasin_TWh','energyReqUIB_TWh' , 'energyReqIndus_TWh', 'energyConsumPak_TWh')
%table(yrsWouter',futpopUIB/1e6,futpopIndus/1e6)
disp("Loaded data")

%% GOOD: Loop heatmap MF and FF in one for sustainable potential
npfs=25;
nqfs=25;

%% GOOD: CLEAN Loop heatmap total generation MF and FF separate from historical pf under historical qf
fns = ["mean" "prctile10" "prctile90"]; %fieldnames(LTannual_basin);
fnnames=["Mean" "10^{th} percentile" "90^{th} percentile"];
tfidx=[1:13; 1 14:25];
myr=1;
minmaxmat_total= table(histrcpcornernames(1:13)','VariableNames',{'pfnames'});
figure
for fi=1:length(fns)
    for futtf=1:2
        datatmp=squeeze(HP_TS{seltechsust}.LTannual_basin.(fns(fi)))/1000;
        datain=datatmp(tfidx(futtf,:),tfidx(futtf,:));
        subplot(3,3,(futtf-1)*3+fi)
        imagesc(datain)
        xline([1.5, 5.5,9.5])%,'LineWidth',1.5)
        yline([1.5, 5.5,9.5])%,'LineWidth',1.5)

        set(gca,'xtick',1:npfs,'xticklabel',histshortcornernames,'ytick',1:npfs,...
            'yticklabel',histshortcornernames,'ticklabelInterpreter','none')%,'xticklabelrotation',90)
        set(gca,'YDir','normal')
        title(fnnames(fi))

        %        colormap(flipud(parula)) %- liked parula coz its brighter
        % fix xy labels and ticks
        if fi==1
            addAlphaLabel(futtf,"outside")
            ylabel([tf_full{futtf};"Portfolios"])
            text(-2*ones(1,3),[3.5:4:12.5],strrep(rcpnames," ",""),'HorizontalAlignment','center','fontweight','bold','FontSize',11, 'Rotation',90);

        elseif fi==2 && futtf==2
            xlabel(["Future discharge scenarios"])%; tf_full{futtf}])
        end

        if fi>1
            yticklabels([])
        end
        % keep xlabel and add rcpnames only in bottom plots below corners

        if futtf==1
            xticklabels([])
        else
            text([3.5:4:12.5],-2*ones(1,3),strrep(rcpnames," ",""),'HorizontalAlignment','center','fontweight','bold','FontSize',11);
        end

        axis square
        caxis([70 140]) % for sustpot

        % Find r,c for best pf for each Q
        [~,R_forPmax]=max(datain,[],1);  % returns max for each column in input matrix
        [~,C_forQmax]=max(datain,[],2);
        [~,R_forPmin]=min(datain,[],1);  % returns max for each column in input matrix
        [~,C_forQmin]=min(datain,[],2);

        % save ridx for best and worst pf for full matrix size
        bestpf_idx(futtf,:,fi)=tfidx(futtf,R_forPmax);
        worstpf_idx(futtf,:,fi)=tfidx(futtf,R_forPmin);

        hold all
        l(1)=scatter(1:13,R_forPmax,"rx",'DisplayName',"Best portfolio for each discharge");
        l(2)=scatter(1:13,R_forPmin,20,"kx",'DisplayName',"Worst portfolio for each discharge");


    end
    % get min max for each pf
    minmaxmat_total.(strcat(tf_full{futtf},"_",fns(fi)))= [min(datain,[],2), max(datain,[],2)];

end
legend(l,'Location','northoutside')%,'NumColumns',2)

colormap(viridis)
% add colorbar at the bottom
cb= mycbar("Annual energy generation in TWh/year", 'southoutside');
cb.Position=[0.09 0.06 0.9 0.015];    % shift it to bottom of fig

%sgtitle(sprintf("annual energy for hist+futpf under hist+futQs for %s potential",currentpotname) );

%% FINAL: TS of energy demand vs supply climate impact vs pf impact

c_4corner=brighten(hsv(4),-0.8); %vga(6); %
figure
histidx=1;
for pot=1:2
    for ssel=1:3
        %        subplot(3,1,ssel)
        subplottight(3,3,ssel)
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


        %         % add energy demand metrics
        %         if pot==2
        %             l1(li+1)=plot(yrsWouter,energyReqIndus_TWh(:,ssel),"--b",'LineWidth',1.5);
        %             l1(li+2)=plot(yrsWouter,energyReqUIB_TWh(:,ssel),":b",'LineWidth',1.5);
        %             l1(li+3)=plot(energyConsumPak_TWh{:,1},energyConsumPak_TWh{:,2},"Color",mybrown,'LineWidth',1); % Historical gen
        %             l1(li+4)=plot(energyConsumPak_TWh{:,1},energyConsumPak_TWh{:,ssel+2},":","Color",mybrown,'LineWidth',1.2); % Future gen projections
        %         end

        applymyplotformat(rcpnames{ssel})

        % fix xy labels and ticks
        xlim([2000, 2095])
        ylim([10 500])
        %xticklabels([])
        if ssel>1
            yticklabels([])
        elseif ssel==1
            ylabel("Hydropower generation (TWh/yr)")
            addAlphaLabel(1,"outside")
        end

        % midfar future labels
        xline([2036, 2066],'LineStyle',':')
        topy=ylim;
        text([2050 2080],[1 1]*floor(topy(2)-50),tname,'HorizontalAlignment','center','FontWeight','bold')
    end
end
%lgd=legend(l1,["Historical", cornernames],'Orientation','horizontal');
%lgd.Title.String="Porfolios";

% FINAL: TS of number of ppl that can be supported by annual TS for the pfs in their own qfs
% hist_annual_ts is 25pfs x 40 yrs, fut is 25pfs x 30 yrs x 25 qs
for pot=1:2
    for ssel=1:3
        %        subplot(3,1,ssel)
        subplottight(3,3,3+ssel)
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
            ylabel("Population (Millions)")
            addAlphaLabel(2,"outside")

        elseif ssel==2
            xlabel("Year")
        end
        if ssel>1
            yticklabels([])
        end
        l12(li+1)=plot(yrsWouter,futpopIndus(:,ssel)/10^6,":b",'LineWidth',1.5);
        l12(li+2)=plot(yrsWouter,futpopUIB(:,ssel)/10^6,":m",'LineWidth',1.5);
    end
end

lgd=legend(l12,["Historical", cornernames, ...
    "Indus ", ...
    "Upper Indus"],'NumColumns',5,'Orientation','horizontal');
lgd.Title.String="Porfolios   Projected Population: ";

%%
disp("***************************************EOF***************************************")

