%% Evaluate actual energy generation at optimal HP projects under historical and future climate
% Converting Q from cell2mat is memory intensive. has tobe done in the HPC.
% Load and prepare portfolios for hist and future scenarios. Loads all not
% nan projects for both tech and sust scenarios
% Check that portfolios are loaded correctly by recalculating potential and
% ensuring it matches the outputs of Hydrus
% Loop through each pot type, load qf, then load each pf and run it for the qf
% Run each portfolio through hist and future Q time series (QmonTS is in m3/mon! so needed to change the units properly)
% Evaluate actual LF and recalculate energy generated
% Plot all pfs in each Q.
% tech and sust pfs saved in separate files due to large size

% Converting Q500m from cell2mat is memory intensive. Code has to be run in the HPC.

clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

loadPortfolios=01; % load saved portfolios
ofmatPortfolios=fullfile(rootof,'FutScenarios_Portfolios02.mat');
ofmatHPgenTS=fullfile(rootof,'FutScenarios_HPgen_allTS03');

% Create tech/fin names for use w Tech and Sust portfolios  
mypotname=[newpots4([1,2]); newpots4([3,4])];
logf = fullfile(rootoffut, sprintf('Run HPsimulation03 on %s.log',string(datetime('now','Format','dd-MMM-yyyy'))));
diary(logf)
tic

%% Extract portfolio or load them from saved ones from previous compile
if loadPortfolios
    load(ofmatPortfolios,'portfolios','perfStats','cl_Qnames','cl_pfnames')
    disp("Loaded HP portfolios and filenames for PF and Q_TS files")
else

    %% Get file names for hist/fut portfolios and hist/fut Q_TS
    % HP portfolios
    load(fullfile(rootof,'MainFutScenarios.mat'),'newfpaths')
    newfpaths=fullfile(rootof,'Full', 'ASIA_Historical_Runs',extractAfter(newfpaths,'ASIA\'));
    portfolio_fpaths=[fullfile(rootof,'Full', 'ASIA_Historical_Runs',["R103_Energy_Full_Mixed_Tech_Fin"; "R103_Energy_Full_Mixed_Sust_RiskAverse"]);
        newfpaths];
    cl_pfnames=strrep(extractAfter(portfolio_fpaths,"R103_"),"Full_Mixed_","");
    cl_pfnames(1:2)=["Historical-Tech" "Historical-Sust"];
    nportfolios=length(cl_pfnames);

    % Q_TSs
    Q_TS_fnames=["MLAB500m_FutQ_TS_m3day_hist.mat";
        strcat("MLAB500m_FutQ_TS_m3day_",[strcat(modorder.matlabname,"_monTS_",tframe{1}); strcat(modorder.matlabname,"_monTS_",tframe{2})],".mat")];
    cl_Qnames=strcat("Q_",strrep(extractAfter(Q_TS_fnames,"_m3day_"),".mat",""));
    cl_Qnames(1)="Q_Historical";
    figure; hold on
    for f=1:nportfolios
        mydata=load(fullfile(portfolio_fpaths{f},sprintf('Basin%d_output_do.mat', nbasin)),...
            'PID','ross','coss','COEAlls','PnetAlls','SysIDs','ross','coss','latss','lonss','enablesmallPcost',...
            'QDesignPinletMinend','QDesignMeanPinletMinend','QDesignLFPinletMinend', 'HeadPminend', 'PPMinend','PPnetend','PPtheoryend','OpthfPminend', 'OptDminend','nPipePMinend','dfPLminEnd',...
            'Q_RD_design','OptDH','RDP','RDPnet','Q_RD_design_LF');

        %From model code, I know that DP data is above RP data COEAll = [horzcat(PCOEend{:})';COETotRD]; PnetAll = [horzcat(PPnetend{:})';RDPnet];
        selDP=1:length(horzcat(mydata.QDesignPinletMinend{:})');
        selRP=selDP(end)+1:length(mydata.COEAlls);
        portfoliosize=size(mydata.COEAlls);

        % data already joined for DP and RP
        p_r=mydata.ross;  %For DP this tracks inlets. For RP this tracks outlets.
        p_c=mydata.coss;
        COEalls=mydata.COEAlls;
        SysIDs=mydata.SysIDs;
        %Pnetalls_GWhs=mydata.PnetAlls;
        Pnetalls_GWhs=[horzcat(mydata.PPnetend{:})'; mydata.RDPnet];  % these are in GWh.
        %CHECKED that this is same as Pnetalls but decided to use this as
        %it has energy even for not selected projects

        % data not yet joined for DP and RP
        Qdesigns=[horzcat(mydata.QDesignPinletMinend{:})'; mydata.Q_RD_design(:)];
        Qdesigns_mean=[horzcat(mydata.QDesignMeanPinletMinend{:})'; mydata.Q_RD_design(:)]; % I do not save Qmeans - LT annualmonavgs for RP so keeping Qdesigns here
        Hgross=[horzcat(mydata.HeadPminend{:})'; mydata.OptDH(:)]; % elevation diff for DP and dam height for RP
        MWs=[horzcat(mydata.PPMinend{:})';   mydata.RDP(:)]/1e6; % W to MW
        LF=[horzcat(mydata.QDesignLFPinletMinend{:})'; mydata.Q_RD_design_LF(:)];
        hfriction=[horzcat(mydata.OpthfPminend{:})'; zeros(size(mydata.RDP(:)))].*Hgross; % only valid for DP, so added 0s to RP. hf is saved in terms of % of Hgross so need to get it back

        % Eval derived params
        % Eval Hnet, mostly for DP than RP
        Hnet=Hgross-hfriction;
        run('config_costmodels.m')
        % assume everything is large first then reset using
        eta_generation = [ones(size(selDP))'*eta_generation_DP_large; ones(size(selRP))'*eta_generation_RD];
        eta_transmission = [ones(size(selDP))'*eta_transmission_DP_large; ones(size(selRP))'*eta_transmission_RD];
        plant_size=2*ones(size(eta_transmission)); % 2=Large, 1=Small

        enablesmallPcost=mydata.enablesmallPcost;

        %% For DP: Classify into small projects for DP based on theory pot
        %PTgross_MW=nan(portfoliosize);
        PTgross_MW = rho*g*Qdesigns_mean(selDP).*Hgross(selDP)*1e-6;    % Theoretical gross potential in MW
        % get eta for small projects
        if enablesmallPcost
            selsmallDP=find(PTgross_MW <50);
            plant_size(selsmallDP)=2;
            eta_generation(selsmallDP) = eta_generation_DP_small;
            eta_transmission(selsmallDP)= eta_transmission_DP_small;
        end

        %% Calculate energy for RP and DP
        Pdesign_GWh=nan(portfoliosize);
        myTheoryPnet_GWh = rho*g*Qdesigns_mean.*Hnet*8760*1e-9; % Theoretical net potential
        myPror_MW = eta_generation.*eta_transmission.*rho*g.*Hnet.*Qdesigns*1e-6;            % Capacity based on calculated Q (MW)
        Pdesign_GWh = (myPror_MW*1e-3*8760).*LF;                   % Yearly production with LF (kWh)

        % Plot difference between my and hydrus calculations
        plot(Pdesign_GWh, Pnetalls_GWhs,'*-','DisplayName',cl_pfnames{f});hold on
        plot(Pdesign_GWh(selRP), Pnetalls_GWhs(selRP),'o')
        plot(Pdesign_GWh(selDP), Pnetalls_GWhs(selDP),'s');grid on
        % legend(clnames{f},"OnlyRP","OnlyDP",'Interpreter','none')
        xlabel("My power gen")
        ylabel("Hydrus power gen")

        %         figure;plot(Pdesign_GWh,'+');hold on
        %         plot(Pnetalls_GWhs,'o')
        %         disp('Evaluated both energy')

        %% Cleanup portfolio
        selvalid=~isnan(COEalls);
        portfolios{f}=table(mydata.PID(selvalid),p_r(selvalid), p_c(selvalid),SysIDs(selvalid),plant_size(selvalid),COEalls(selvalid), ...
            Pnetalls_GWhs(selvalid),MWs(selvalid),Pdesign_GWh(selvalid),Qdesigns_mean(selvalid), ...
            Qdesigns(selvalid),LF(selvalid),Hnet(selvalid),eta_generation(selvalid),eta_transmission(selvalid), ...
            'VariableNames',{'pid','r_idx','c_idx','SysID','Plant_size','COEall','Pnetall_GWh','MW','Pdesign_GWh','Qdesign_mean','Qdesign','LFdesign','Hnet','eta_generation', 'eta_transmission'});

        % portfolios{f}.tmpPactual_Wh = portfolios{f}.eta_generation.*portfolios{f}.eta_transmission.*rho*g.*portfolios{f}.Hnet; % this needs to be multiplied by Qt*1e-9 for HPgen = eta*rho*g*H*(Qdesign*LF*t) *1e-9 to get energy in GWh/yr

        perfStats(f,:)=[nansum(Pdesign_GWh-Pnetalls_GWhs) nansum(Pdesign_GWh(selRP)-Pnetalls_GWhs(selRP)) nansum(Pdesign_GWh(selDP)-Pnetalls_GWhs(selDP))];
    end
    legend("Interpreter","none")

    %% Archive portfolio
    save(ofmatPortfolios,'portfolios','perfStats','cl_Qnames','cl_pfnames')
    disp('Saved portfolios and filenames')
end

%% Create timeseries for days in month separately for hist and future
hist_daysinmon_ts=eomday(year(histmon_ts),month(histmon_ts));
fut1_daysinmon_ts=eomday(year(futmon_ts1),month(futmon_ts1));
fut2_daysinmon_ts=eomday(year(futmon_ts1),month(futmon_ts1));

%% Get HP generation for each discharge TS and portfolio combination
%QF files are ordered hist followed by ssps for near and the far future
%PF files are ordered w hist Tech/Sust followed near future ssp245 Tech/Sust and then far future diff ssps
%create indices to only run tech first then sust - 25 tech, 25 sust
techidx=1:2:50;
sustidx=2:2:50;
pfgroups= [techidx; sustidx];
%reorder pf names to have techidx first
pfnames=cl_pfnames(pfgroups)';

eflow=30/100;           %Percentage of water diverted around the dam to ensure natural ecological flow of the river
nqfs=length(cl_Qnames);

% Loop through each pot type, load qf, then load each pf and run it for the qf
% need to do pot before Q so each pot can be saved one at a time though
% reloading Q for each pot is a bit slow.
for pot= 2:-1:1
    % Loop through each discharge time series and eval HP for all portfolios
    % under each Q
    for qf=1:nqfs
        % Load Q in m3day and convert to m3s
        if qf==1
            QTSmatpath=fullfile(rootf,"data","UI","data","Q_monTS_m3month","MLAB500m_FutQ_TS_m3day_hist.mat");
        else
            QTSmatpath=fullfile(rootf,"data","UI","data","Q_monTS_m3month",strcat("MLAB500m_FutQ_TS_m3day_",extractAfter(cl_Qnames{qf},"Q_"),".mat"));
        end
        load(QTSmatpath,'Q500m_m3day') % <<< THIS DATA IS ACTUALLY in m3/month

        % Apply eflow only for sust potential, note that in small HPs eflow
        % is not applied in HyPE so should technically be done here too but
        % i dont have a tacking of small vs large DP
        if pot==2
            Q_TS_m3month=cat(3,Q500m_m3day{:})*( 1-eflow); % convert cell to matrix
        else
            Q_TS_m3month=cat(3,Q500m_m3day{:}); % convert cell to matrix
        end
        clear Q500m_m3day
        fprintf("Calculating HPgeneration TS for %s pot under %s #%d out of %d Qs\n",pottypes2{pot},cl_Qnames{qf},qf,nqfs)

        % sel right tf for fut - mid or far
        if qf<=13 %MF
            sel_fut_daysinmon_ts=fut1_daysinmon_ts;
        else %FF
            sel_fut_daysinmon_ts=fut2_daysinmon_ts;
        end

        % Loop through each portfolio
        pf=0;
        for pfi=pfgroups(pot,:)
            pf=pf+1;

            % Loop through each project
            nprojects=height(portfolios{pfi});
            for pt=1:nprojects
                %  pt=7330 and 100
                r0=portfolios{pfi}.r_idx(pt);
                c0=portfolios{pfi}.c_idx(pt);

                if qf==1 %histQ has more timesteps than futQ so cant save it in the same structure so save it separately
                    % Extract Q TS for all plant locations - inlet for DP and dam for HP
                    actualQ_m3s=squeeze(Q_TS_m3month(r0,c0,:))./hist_daysinmon_ts/(24*60*60); % Convert m3/month to m3/s
                    actualLF=min(1,actualQ_m3s/portfolios{pfi}.Qdesign(pt));

                    % Energy = Power in MW * LF * 24hs * ndays in month where MW already = eta*rho*g*H*Qdesign
                    PactualTS_GWh_4histQ{pf}(pt,:,qf) = portfolios{pfi}.MW(pt).*actualLF*24.*hist_daysinmon_ts*1e-3 ;  % HPgen saved in a 3d mat w projects x Qtimeseries x Q_scenarios
                else

                    % Extract Q TS for all plant locations - inlet for DP and dam for HP
                    actualQ_m3s=squeeze(Q_TS_m3month(r0,c0,:))./sel_fut_daysinmon_ts/(24*60*60); % Convert m3/month to m3/s
                    actualLF=min(1,actualQ_m3s/portfolios{pfi}.Qdesign(pt));

                    % Get monthly timeseries of energy generated
                    % Energy = Power in MW * LF * 24hs * ndays in month where MW already = eta*rho*g*H*Qdesign
                    PactualTS_GWh_4futQ{pf}(pt,:,qf) = portfolios{pfi}.MW(pt).*actualLF*24.*sel_fut_daysinmon_ts*1e-3 ;  % HPgen saved in a 3d mat w projects x Qtimeseries x Q_scenarios
                end
            end
            % sel only financial pot
            selprjs=find(portfolios{pfi}.COEall<=costlim)';
            if qf==1
                PactualTS_GWh_4histQ_fin{pf}(:,:,qf) = PactualTS_GWh_4histQ{pf}(selprjs,:,qf);
            else
                PactualTS_GWh_4futQ_fin{pf}(:,:,qf) = PactualTS_GWh_4futQ{pf}(selprjs,:,qf);
            end

            % Get basin level total
            if qf==1
                tmp_Pactual_monLT =reshape(sum(PactualTS_GWh_4histQ{pf}(:,:,qf))/1000,12,[]);
                Pactual_annualTS_4histQ(pf,:,qf)=sum(tmp_Pactual_monLT);
            else
                tmp_Pactual_monLT =reshape(sum(PactualTS_GWh_4futQ{pf}(:,:,qf))/1000,12,[]);
                Pactual_annualTS_4futQ(pf,:,qf)=sum(tmp_Pactual_monLT);
            end
        end
    end
   

    %% Plot annual TS of all pfs under histQ
    figure
    plot(Pactual_annualTS_4histQ')
    legend(cl_pfnames(techidx),'Location','eastoutside','Interpreter','none')
    title(mypotname(pot,1))
    hold on
    plot(Pactual_annualTS_4histQ(1,:),'k','LineWidth',2,"DisplayName","Historical")

    %% Archive HP generation - all pf and only fin
    save(strcat(ofmatHPgenTS,"_",mypotname(pot,1),".mat"),...
        'PactualTS_GWh_4futQ','PactualTS_GWh_4histQ','cl_Qnames','pfnames','-v7.3')
    save(strcat(ofmatHPgenTS,"_",mypotname(pot,2),".mat"),...
        'PactualTS_GWh_4futQ_fin','PactualTS_GWh_4histQ_fin','cl_Qnames','pfnames','-v7.3')
    fprintf("Saved HP generation TS stats for %s potential\n", mypotname(pot,1:2))
    clear Pactual*
end
toc
diary off