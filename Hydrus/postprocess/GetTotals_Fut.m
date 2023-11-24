% Post-process to get total energy and # of projects for FUTURE scenario
% runs. Loops through folders based on folder name pattern specified
% option to write separate excel files for CAS
% Save basin level total Tech, Fin, Econ potential for all scenarios. For
% Creates generic bar plots for all scenarios
% Compares changes in DP and RP potentials
% Creates cleaned up bar plots for:
        % Horizontal bar plot for basin level pot tots
        % Vertical bar plot for basin level pot tots and num projs
        % change in climate and hydrology vs change in HP

clc
clearvars
close all
% Add folder path
addpath('Hydrus\postprocess')
run('myVarNames_Fut.m')
% varname
loadsaved=10;
printsummary=0;
plotoldfigs=0;
writeSummary2file=0;
write4CAS=0;
% output fnames
ofname='Compile_FutScenario_Totals00';
xlsfile_summary=fullfile(rootof,[ofname '.xlsx']);
matfile_summary=fullfile(rootof,[ofname '.mat']);
xlspath=fullfile(rootf,'CAS_19July2022'); if ~isfolder(xlspath); mkdir(xlspath); end
xlsVarnames={'pid', 'Latss', 'Lonss', 'SysIDs', 'COEAlls_[USD2010/kWh]', 'PnetAlls_[GWh/yr]'};

addpath(genpath(fullfile(rootf,"..","Hydrus")))

%% Load future tech, fin, sust data - save or load here
if loadsaved
    load(matfile_summary,'tot_allscen','potlist','histHP')
    disp("All scen tots loaded")
else
    %% Create correct order of models (48 = 4cornersx 3SSP x 2TF x 2pottype)
    modorder=readtable(fullfile(rootoffut,"FutScen_Tracker.xlsx"),Sheet="GCMdetails", Range="A1:F13");
    suffix=strrep(tframe,'-','_');
    suffix=[strcat("_",suffix(1),["_Tech_Fin" "_Sust_RiskAverse"]) strcat("_",suffix(2),["_Tech_Fin" "_Sust_RiskAverse"])]';
    myforder=[];

    % Create neworder to have for each sspmodel in one place
    for f=1:height(modorder)
        myforder=[myforder ;strcat(modorder.matlabname(f),suffix)];
    end

    % Get filepaths
    [fpaths, fnames]=path2fldrfiles(fullfile(rootf),strcat(r_prefix,"*"));
    % match order to new
    newforder=zeros(1,length(myforder));
    fnames_cl=strrep(strrep(strrep(fnames,"FutR103_",""),"_Full_Mixed",""),"_LTavgs","");
    for f=1:length(myforder)
        newforder(f)=find(fnames_cl==myforder(f));
    end
    newfnames=fnames_cl(newforder');
    newfpaths=fpaths(newforder);
    disp("Reorderd future scenarios")

    %% Loop through and compile FULL + REMAIN 3 policy types
    idx=1; % idx to track tech, fin, sust scenarios
    for fidx=1:length(newfpaths)
        runname=newfnames{fidx};
        load(fullfile(newfpaths{fidx},sprintf('COEPOT_b%d_do.mat', nbasin)))

        % Eval tech and econ data
        seltech= COEAlls>=0;
        selecon= COEAlls<=costlim;
        tech_outdata_ss=table(PID(seltech),latss(seltech),lonss(seltech),SysIDs(seltech),COEAlls(seltech),PnetAlls(seltech),...); % exported to excel
            'VariableNames',xlsVarnames);
        fin_outdata_ss=table(PID(selecon),latss(selecon),lonss(selecon),SysIDs(selecon),COEAlls(selecon),PnetAlls(selecon),...); % exported to excel
            'VariableNames',xlsVarnames);
        if printsummary
            fprintf('\n\nFor %s\n', runname)
            fprintf('%s Technical Potential: %0.2f GWh\n',scenario,nansum(PnetAlls))
            fprintf('# of selected projects: %0.0f out of %0.0f\n',sum(~isnan(PnetAlls)),length(PnetAlls))
            fprintf('# of selected River type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==2),sum(PnetAlls(SysIDs==2)))
            fprintf('# of selected Diversion type projects: %0.0f with %0.2f GWh\n\n',sum(SysIDs==1),sum(PnetAlls(SysIDs==1)))
            fprintf('%s Financial Potential: %0.2f GWh\n',scenario, nansum(PnetAlls(COEAlls<=costlim)))
            fprintf('# of selected projects: %0.0f out of %0.0f\n',sum(~isnan(PnetAlls(COEAlls<=costlim))),length(PnetAlls))
            fprintf('# of selected River type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==2 & COEAlls<=costlim),sum(PnetAlls(SysIDs==2 & COEAlls<=costlim)))
            fprintf('# of selected Diversion type projects: %0.0f with %0.2f GWh\n\n',sum(SysIDs==1 & COEAlls<=costlim),sum(PnetAlls(SysIDs==1 & COEAlls<=costlim)))
        end

        % Evaluate % of projects below <=costlim in terms of energy
        prctfin(idx,1) = sum(~isnan(PnetAlls(COEAlls<=costlim)))/sum(~isnan(PnetAlls))*100;
        prctfin(idx+1,1) = 100; % for financial potential part
        %Archive
        compile_PNetalls(idx,:)= [sum(PnetAlls(SysIDs==2)) sum(PnetAlls(SysIDs==1)) sum(SysIDs==2) sum(SysIDs==1)]; %Tech - RP/DP GWh/yr and counts
        compile_PNetalls(idx+1,:)= [sum(PnetAlls(SysIDs==2 & COEAlls<=costlim)) sum(PnetAlls(SysIDs==1 & COEAlls<=costlim)) sum(SysIDs==2 & COEAlls<=costlim) sum(SysIDs==1 & COEAlls<=costlim)]; %Fin - RP/DP GWh/yr and counts
        compile_runname{idx,1}=[runname '_Tech'];
        compile_runname{idx+1,1}= [runname '_Fin'];
        idx=idx+2;

        %% Write Tech/Fin to excel
        if write4CAS
            oxlsfile=fullfile(xlspath,[runname '.xlsx']);
            writetable(tech_outdata_ss,oxlsfile,'Sheet','Tech_RiskAverse')
            writetable(fin_outdata_ss,oxlsfile,'Sheet','Fin_RiskAverse')
        end
    end
    disp("Loaded all future scenarios")
    %% Compile all scenarios and their details
    tot_allscen=table(compile_runname, compile_PNetalls(:,1)/1000, compile_PNetalls(:,2)/1000, sum(compile_PNetalls(:,1:2),2)/1000,...
        compile_PNetalls(:,3), compile_PNetalls(:,4), sum(compile_PNetalls(:,3:4),2), prctfin,...
        'VariableNames',{'Runname','TWh/yr of River power','TWh/yr of Diversion canal','Total TWh/yr','Number of River power',...
        'Number of Diversion canal','Total number','Prct nplants costing <=0.10'});

    % Get run details
    fdetails=split(compile_runname,'_');
    cornerlist=categorical(repmat(repelem(cornernames,8),1,3)',cornernames,'Ordinal',1);
    potlist=categorical(join(fdetails(:,[5,7]),"-"),newpots4,'Ordinal',1);
    scen_details = table(categorical(fdetails(:,1),unique(fdetails),'Ordinal',1),cornerlist,fdetails(:,2),join(fdetails(:,3:4),"-"),potlist, ...
        'VariableNames',{'ssp','corner','gcm','yrs','pottype'});
    alldata=[tot_allscen scen_details];
    % reorder so scens are sorted by
    tot_allscen=sortrows(alldata,{'yrs','ssp','pottype','corner','Total TWh/yr'},{'ascend','ascend','ascend','ascend','descend'});
    disp("Sorted and compiled all data")

    %% Add empty rows to force space between three energy scens
    tot_allscen{end+1,1}={' '};    tot_allscen{end,2:8}=nan;

    %% Load and save only mixed scenario of historical tech, fin, sust potential
tmp=load(fullfile(rootof,'Compile_Scenario_Totals03.mat'), 'tot_allscen');
selhist=7:9; % Mixed scen tech, fin and sust-risk averse
histHP=tmp.tot_allscen(selhist,:);
clear tmp
% manually add the Sust-tech part
histHP(4,:)= {'Energy-Full-Mixed-Sust-Tech' 61.755 46.936 108.691 220 4033 4253};
% reorder so sust tech is moved up
histHP=histHP([1,2,4,3],:);
histHP.("Prct nplants costing <=0.10") = [histHP.("Total number")(2)/histHP.("Total number")(1) 1 histHP.("Total number")(4)/histHP.("Total number")(3) 1]'*100;
histHP.pottype=newpots4';
disp("Loaded historical mixed scen data")
end

%% Load hist and future theoretical potential
load(fullfile(rootoffut,'FutTheoreticalPot_subbasintotals.mat'), 'subPot_GWh','subPot_prctchange',"ssp_gcm_tf_names")

%% Evaluate change in tech/fin/sust potential
change_allscen=tot_allscen; %change in magnitude
changeprct_allscen=tot_allscen; %change in terms of % of historical value

tot_allscen.("Historical TWh/yr")=nan(height(tot_allscen),1);
tot_allscen.("Historical Num Proj")=nan(height(tot_allscen),1);

for p=histHP.pottype'
    selp=change_allscen.pottype==p;
    change_allscen(selp,2:8)= array2table(change_allscen{selp,2:8}-histHP{histHP.pottype==p,2:8});

    % Add hist data to tot_allscen
    tot_allscen.("Historical TWh/yr")(selp)=histHP{histHP.pottype==p,4};
    tot_allscen.("Historical Num Proj")(selp)=histHP{histHP.pottype==p,7};
    tot_allscen.("Historical Prct nplants costing <=0.10")(selp)=histHP{histHP.pottype==p,8};
end

%evaluate change in terms of % of historical value
changeprct_allscen{:,2:4}=change_allscen{:, 2:4}./tot_allscen{:,"Historical TWh/yr"}*100;
changeprct_allscen{:,5:7}=change_allscen{:, 5:7}./tot_allscen{:,"Historical Num Proj"}*100;
changeprct_allscen{:,8}=change_allscen{:, "Prct nplants costing <=0.10"}./tot_allscen{:,"Historical Prct nplants costing <=0.10"}*100;
disp("Evaluated change in HP")

%% Compare RP vs DP
% skip sust-tech case
selrows=change_allscen.pottype~="Sust-Tech";
%selrows=1:height(change_allscen); %

%RP-DP
diffmagRPvsDP=change_allscen{selrows,2}-change_allscen{selrows,3};
diffnRPvsDP=change_allscen{selrows,5}-change_allscen{selrows,6};


fprintf("Magnitude Change in DP energy > RP energy in %0.3f %% of scenarios across tech/fin/sust pottypes\n",sum(change_allscen{selrows,3}>change_allscen{selrows,2})/height(change_allscen)*100)
fprintf("Magnitude Change in DP number > RP number in %0.3f %% of scenarios across tech/fin/sust pottypes\n\n",sum(change_allscen{selrows,6}>change_allscen{selrows,5})/height(change_allscen)*100)

% % change does not make sense to calculate actually becuase this i %
% change in full hist pot, not % change in DP or RP specifically
%fprintf("%% Change in DP energy > RP energy in %0.3f %% of scenarios across tech/fin/sust pottypes\n",sum(changeprct_allscen{selrows,3}>changeprct_allscen{selrows,2})/height(changeprct_allscen)*100)
%fprintf("%% Change in DP number > RP number in %0.3f %% of scenarios across tech/fin/sust pottypes\n",sum(changeprct_allscen{selrows,6}>changeprct_allscen{selrows,5})/height(changeprct_allscen)*100)

% Eval change as % of RP/DP
for pt=newpots4([1,2,4])
    selp=find(change_allscen.pottype==pt);
    hist_selrow=find(histHP.pottype==pt);

    change_allscen.dRPTWh_asprctRP(selp)=change_allscen{selp,2}/histHP{hist_selrow,2}*100;
    change_allscen.dRPnprj_asprctRP(selp)=change_allscen{selp,5}/histHP{hist_selrow,5}*100;

    change_allscen.dDPTWh_asprctDP(selp)=change_allscen{selp,3}/histHP{hist_selrow,3}*100;
    change_allscen.dDPnprj_asprctDP(selp)=change_allscen{selp,6}/histHP{hist_selrow,6}*100;
end

% show min/max of RP and DP change in total vs change in
tmp_min=min(change_allscen{selrows,[2:6, 17, 19, 18, 20]});
tmp_max=max(change_allscen{selrows,[2:6, 17, 19, 18, 20]});

RPDPchange=array2table([tmp_min(1) tmp_max(1) tmp_min(4) tmp_max(4) tmp_min(6) tmp_max(6) tmp_min(8) tmp_max(8)
                        tmp_min(2) tmp_max(2) tmp_min(5) tmp_max(5) tmp_min(7) tmp_max(7) tmp_min(9) tmp_max(9)], ...
    'VariableNames', {'min TWhperyr', 'max TWhperyr', 'min Dnumprj', 'max Dnumprj', '% min TWhperyr', '% max TWhperyr', '% min Dnumprj', '% max Dnumprj'},'RowNames',planttypes(:))


disp("Compared changes in RP and DP HP")

%% Plot DP RP changes
figure;
subplot(4,1,1)
bar(change_allscen{selrows,2:3})
applymyplotformat("Change in TWh/yr for RP DP tech/fin/sust potentials across all fut scenarios")
subplot(4,1,2)
bar(change_allscen{selrows,5:6})
applymyplotformat("Change in # of project for RP DP tech/fin/sust potentials across all fut scenarios")
%set(gca,'YScale','log')
colororder([cl_RP;cl_DP])
legend(planttypes)

subplot(4,1,3)
plot(diffmagRPvsDP,"*")
ylabel("\Delta TWh/yr")

hold all
yyaxis right
plot(diffnRPvsDP,"+")
ylabel("\Delta numproj")
applymyplotformat("dRP-dDP")
legend("Change in energy","Change in numprj")

subplot(4,1,4)
plot(change_allscen{selrows,[17, 19]},"*")
ylabel("% \Delta in Energy TWh/yr")

hold all
yyaxis right
plot(change_allscen{selrows,[18, 20]},"o")
ylabel("% \Delta in numproj")
applymyplotformat("% changes in RP vs DP")
legend("Change in energy","Change in numprj")

%% Figs
if plotoldfigs
    %% GOOD: Plot bar for 96 future energy scenarios
    figure
    nscens=3*4*2*4;
    selfull=[1:16 97 17:32 97 33:48 97 49:64 97 65:80 97 81:96];% add idx 97 as cushion
    xx=[2.5:4:15, 19.5:4:31.5 36.5:4:50 53.5:4:67 70.5:4:82.5 87.5:4:100];

    % Plot bar pts for totalpot
    b1=bar(tot_allscen{selfull, 2:3},.95,'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',baralpha,'EdgeAlpha',0);
    hold all

    % Change color of  bar
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};

    end

    % Add labels for year and pot type
    ylabel('Energy in TWh/year','FontWeight','bold')
    xticks([])

    yyaxis right
    scshift=0.1;
    s1(1)=scatter([1:length(selfull)]-scshift,tot_allscen{selfull, 5},30,cl_RP*.6,'d','filled');
    hold on
    s1(2)=scatter([1:length(selfull)]-scshift,tot_allscen{selfull, 6}+tot_allscen{selfull, 5},30, cl_DP*.6,'d','filled'); %manually stack scatter
    ylabel('Number of projects','FontWeight','bold')
    ax = gca;
    ax.YAxis(2).Color = 'k';
    xticks([])

    grid on
    set(gca,'XGrid','off')
    % label time frames
    text(length(selfull)/4,12200,tframe(1),'FontWeight','bold','FontSize',16,'HorizontalAlignment','center')
    text(length(selfull)/4*3,12200,tframe(2),'FontWeight','bold','FontSize',16,'HorizontalAlignment','center')

    %label potential types
    text(xx(1:4),repelem(11800,4),tot_allscen.pottype(1:4:16),'HorizontalAlignment','center','FontAngle','italic')
    % label ssps
    text(xx([2:4:24])+2,-300*ones(1,6),repmat(modelnames(3:5),1,2),'HorizontalAlignment','center','FontWeight','bold','FontAngle','italic')
    % label scenarios
    %set(gca,'XTick',1:102,'xticklabels', tot_allscen.Runname(selfull),'TickLabelInterpreter','none','XTickLabelRotation',90)

    % add separator lines
    xline(51,'k','LineWidth',2) % for tf
    xline([17 34 68 85],'Color',mygraylines,'LineWidth',1.15)  % for ssp
    xline(xx([1:3,5:7,9:11,13:15, 17:19, 21:23])+2,'Color',mygraylines,'LineStyle',':','LineWidth',1.15)

    %i=1; %add fig label
    %text(0.01,300,char(64+i) ,'FontWeight','bold','FontSize',16)

    yyaxis left
    % Plot historical
    l1=plot(tot_allscen.("Historical TWh/yr")(selfull),'_b','MarkerSize',10);

    % label corners
    text(1:5-.5,tot_allscen{1:4,"Total TWh/yr"}+10,strcat("\leftarrow ",cornernames),'Rotation',90,'FontAngle','italic') %repelem(350,4) %

    % add legend
    legend([b1 s1 l1],["RP: Energy", "DP: Energy", ...
        "RP: Number","DP: Number","Historical Energy"], "NumColumns", 3,'Location','northoutside') %,'Orientation','horizontal')%)
    l.Title.String=sprintf("Energy Potential | Number of Projects ");

    %% GOOD: Plot bar for change in fut energy
    figure
    idata=change_allscen;
    % Plot bar pts for totalpot
    b1=bar(idata{selfull, 2:3},.95,'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',baralpha,'EdgeAlpha',0);
    hold all

    % Change color of  bar
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
        b2(k).EdgeColor = mycolors{k}*0.6;
    end

    % Add labels for year and pot type
    ylabel('Change in Energy in TWh/year','FontWeight','bold')
    xticks([])

    % Plot scattter pts for number on right axis

    yyaxis right
    scshift=0.1;
    scatter([1:length(selfull)]-scshift,idata{selfull, 5},30,cl_RP*.6,'d','filled')
    hold on
    scatter([1:length(selfull)]-scshift,idata{selfull, 6}+idata{selfull, 5},30, cl_DP*.6,'d','filled') %manually stack scatter
    ylabel('Change in Number of projects','FontWeight','bold')
    ax = gca;
    ax.YAxis(2).Color = 'k';
    xticks([])
    ylim([-1000 4000])

    grid on
    set(gca,'XGrid','off')
    % label time frames
    text(length(selfull)/4,4100,tframe(1),'FontWeight','bold','FontSize',16,'HorizontalAlignment','center')
    text(length(selfull)/4*3,4100,tframe(2),'FontWeight','bold','FontSize',16,'HorizontalAlignment','center')

    %label potential types
    text(xx(1:4),repelem(3900,4),idata.pottype(1:4:16),'HorizontalAlignment','center','FontAngle','italic')
    % label ssps
    text(xx([2:4:24])+2,-1100*ones(1,6),repmat(modelnames(3:5),1,2),'HorizontalAlignment','center','FontWeight','bold','FontAngle','italic')
    % label scenarios
    %set(gca,'XTick',1:102,'xticklabels', tot_allscen.Runname(selfull),'TickLabelInterpreter','none','XTickLabelRotation',90)

    % add separator lines
    xline(51,'k','LineWidth',2)
    xline([17 34 68 85],'Color',mygraylines,'LineWidth',1.15)
    xline(xx([1:3,5:7,9:11,13:15, 17:19, 21:23])+2,'Color',mygraylines,'LineStyle',':','LineWidth',1.15)

    %i=1; %add fig label
    %text(0.01,300,char(64+i) ,'FontWeight','bold','FontSize',16)

    % label corners
    yyaxis left
    text(1:4,idata{1:4,"Total TWh/yr"}+10,strcat("\leftarrow ",cornernames),'Rotation',90,'FontAngle','italic','VerticalAlignment','middle') %repelem(350,4) %

    % add legend
    legend(["RP: Energy", "DP: Energy", ...
        "RP: Number","DP: Number"], "NumColumns", 2,'Location','northoutside') %,'Orientation','horizontal')%)
    l.Title.String=sprintf("Energy Potential | Number of Projects ");

    %% Simple plot for total % change
    figure
    bar(change_allscen{selfull, "Total TWh/yr"}./tot_allscen{selfull,"Historical TWh/yr"}*100, ...
        .95,'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',baralpha,'EdgeAlpha',0);
    hold all
    scatter([1:length(selfull)],change_allscen{selfull, "Total number"}./tot_allscen{selfull,"Historical Num Proj"}*100,30,cl_RP*.6,'d','filled')

    %% GOOD: Plot bar for change in fut energy as % of historical totals
    figure
    idata=changeprct_allscen;
    % Plot bar pts for totalpot
    b1=bar(idata{selfull, 2:3},.95,'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',baralpha,'EdgeAlpha',0);
    hold all

    % Change color of  bar
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
        b2(k).EdgeColor = mycolors{k}*0.6;
    end

    % Add labels for year and pot type
    ylabel('% Change from historical total','FontWeight','bold')
    xticks([])

    % Plot scattter pts for number on same axis
    scshift=0.1;
    scatter([1:length(selfull)]-scshift,idata{selfull, 5},30,cl_RP*.6,'d','filled')
    hold on
    scatter([1:length(selfull)]-scshift,idata{selfull, 6}+idata{selfull, 5},30, cl_DP*.6,'d','filled') %manually stack scatter
    grid on
    set(gca,'XGrid','off')

    % label time frames
    text(length(selfull)/4,63,tframe(1),'FontWeight','bold','FontSize',16,'HorizontalAlignment','center')
    text(length(selfull)/4*3,63,tframe(2),'FontWeight','bold','FontSize',16,'HorizontalAlignment','center')

    %label potential types
    text(xx(1:4),repelem(58,4),idata.pottype(1:4:16),'HorizontalAlignment','center','FontAngle','italic')
    % label ssps
    text(xx([2:4:24])+2,-11*ones(1,6),repmat(modelnames(3:5),1,2),'HorizontalAlignment','center','FontWeight','bold','FontAngle','italic')
    % label corners
    text([1:4]-.25,idata{1:4,"Total TWh/yr"}+2,strcat("\leftarrow ",cornernames),'Rotation',90,'FontAngle','italic') %repelem(350,4) %
    % label scenarios
    %set(gca,'XTick',1:102,'xticklabels', tot_allscen.Runname(selfull),'TickLabelInterpreter','none','XTickLabelRotation',90)

    % add separator lines
    xline(51,'k','LineWidth',2)
    xline([17 34 68 85],'Color',mygraylines,'LineWidth',1.15)
    xline(xx([1:3,5:7,9:11,13:15, 17:19, 21:23])+2,'Color',mygraylines,'LineStyle',':','LineWidth',1.15)

    %i=1; %add fig label
    %text(0.01,300,char(64+i) ,'FontWeight','bold','FontSize',16)

    % add legend
    legend(["RP: Energy", "DP: Energy", ...
        "RP: Number","DP: Number"], "NumColumns", 2,'Location','northoutside') %,'Orientation','horizontal')%)
    l.Title.String=sprintf("Energy Potential | Number of Projects ");


end
%% GOOD: Horizonatal bar plot only total

figure
subplot(1,4,1)
% plot mid and far terms as separate
barh(reshape(subPot_GWh(9,2:end)/1000,12,2)','EdgeAlpha',0)
xline(subPot_GWh(9,1)/1000,'color','red','LineWidth',1.5)
yline(1.5,'color','black','LineWidth',1.05) % tf splits
set(gca,'XDir','reverse','YAxisLocation','right','YDir','reverse','YTickLabels',[]) %,'XLim',[0 110])
applymyplotformat('Theoretical',c_ssp)

pp=2;
for pt=newpots4([1,2,4])
    selp=find(tot_allscen.pottype==pt);

    % plot total
    subplot(1,4,pp)
    b1=barh(reshape(tot_allscen.("Total TWh/yr")(selp),12,2)','EdgeAlpha',0);
    hold all

    xline(tot_allscen.("Historical TWh/yr")(selp(1)),'color','red','LineWidth',1.5)
    yline(1.5,'color','black','LineWidth',1.05) % tf splits
    xlim([0 440])
    applymyplotformat(pottypes3(pp-1),c_ssp)
    if pp==2
        l=legend([repelem({' ',},24/2-4) cornernames 'RP' 'DP' 'Hist'],'NumColumns',4);
        l.Title.String=strjoin(rcpnames,'      ');
        yticklabels(strcat(tname,": ", tframe))
    else
        yticklabels([])
    end

    pp=pp+1;
end
%% GOOD: Horizonatal bar plot w total and RP/DP
% Add zero cell to create gaps
subPot_GWh(9,26)=0;

figure
subplot(1,4,1)
% plot mid and far terms as separate
barh(diag(subPot_GWh(9,[2:13,26,14:25])/1000'),'stacked','EdgeAlpha',0)
xline(subPot_GWh(9,1)/1000,'color','red','LineWidth',1.5)
yline(13,'color','black','LineWidth',1.05) % tf splits
set(gca,'XDir','reverse','XColor',[0.00,0.45,0.74],...
    'YDir','reverse','YAxisLocation','left','YTick',1:25,'YTickLabels',[]) %,'XLim',[0 110])
applymyplotformat('Theoretical',[c_ssp;0 0 0])
text([1 1]*3000,[6.5 13+6.5],strcat(tname,": ", tframe),'HorizontalAlignment','center','Rotation',90)
xlim([0 2500])

pp=2;
for pt=newpots4([1,2,4])
    selp=find(tot_allscen.pottype==pt);
    % add gap
    selrows=[selp(1:12); 97; selp(13:24)];
    % plot total and add empty bounds for RP, DP
    subplot(1,4,pp)
    b1=barh(diag(tot_allscen.("Total TWh/yr")(selrows)),'stacked','EdgeAlpha',0);
    hold all
    b2=barh(tot_allscen{selrows, 2:3},.7,  'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.5,'LineStyle','-');
    applymyplotformat(pottypes3(pp-1),[c_ssp;0 0 0])
    set(gca,    'YDir','reverse')
    % Change color of  bar
    for k = 1:2
        b2(k).EdgeColor = mycolors{k}*0.6;
    end

    b3=xline(tot_allscen.("Historical TWh/yr")(selp(1)),'color','red','LineWidth',1.5);
    yticks(1:length(selrows))
    xlim([0 440])

    yline(13,'color','black','LineWidth',1.05) % tf splits

    if pp==2
        l=legend([b1(1:12) b2 b3],[repelem({' ',},12-4) cornernames 'RP' 'DP' 'Hist'],'NumColumns',4);
        l.Title.String=strjoin(rcpnames,'      ');
        yticklabels([rcpcornernames '' rcpcornernames])
    else
        yticklabels([])
    end

    pp=pp+1;
end

%% GOOD: Horizonatal bar plot w total and RP/DP and hist bar
% Add zero cell to create gaps
subPot_GWh(9,26)=0;

figure
subplot(1,4,1)
% plot mid and far terms as separate
sel2plot=[1,26,2:13,26,14:25];
barh(diag(subPot_GWh(9,sel2plot)/1000'),'stacked','EdgeAlpha',0)
xline(subPot_GWh(9,1)/1000,'color','black','LineWidth',1.5,'LineStyle',':')
yline([2,15],'color',mygraylines,'LineWidth',1.05) % tf splits
set(gca,'XDir','reverse','XColor',[0.00,0.45,0.74],...
    'YAxisLocation','left','YTick',1:25,'YTickLabels',[]) %,'XLim',[0 110]),'YDir','reverse',
applymyplotformat('Theoretical',c_ssp26(sel2plot,:))
text([1 1]*3000,[8.5 15+6.5],strcat(tname,": ", tframe),'HorizontalAlignment','center','Rotation',90)
xlim([0 2500])

pp=2;
for pt=newpots4([1,2,4])
    selp=find(tot_allscen.pottype==pt);
    % add gap
    selrows=[97; selp(1:12); 97; selp(13:24)];
    % plot totals
    subplot(1,4,pp)
    idata=[histHP.("Total TWh/yr")(histHP.pottype==pt);   tot_allscen.("Total TWh/yr")(selrows)]; % add historical
    b1=barh(diag(idata),'stacked','EdgeAlpha',0);
    hold all
    % add boxes for RP/DP
    idata=[histHP{(histHP.pottype==pt), 2:3};   tot_allscen{selrows, 2:3}]; % add historical
    b2=barh(idata,.7,  'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.5,'LineStyle','-');
    applymyplotformat(pottypes3(pp-1),c_ssp26(sel2plot,:))
    %set(gca,    'YDir','reverse')
    % Change color of  bar
    for k = 1:2
        b2(k).EdgeColor = mycolors{k}*0.6;
    end

    b3=xline(tot_allscen.("Historical TWh/yr")(selp(1)),'color','black','LineWidth',1.5,'LineStyle',':');
    yticks(1:(length(selrows)+1))
    xlim([0 440])

    yline([2,15],'color',mygraylines,'LineWidth',1.05) % tf splits

    if pp==2
        l=legend([b1(3:14) b2 b1(1) b3],[repelem({'',},12-4) cornernames 'RP' 'DP' 'Historical'],'NumColumns',4, 'Location','northoutside');
        l.Title.String=strjoin(rcpnames,'      ');
        yticklabels([histrcpcornernames(1) '' rcpcornernames '' rcpcornernames])
    else
        yticklabels([])
    end

    pp=pp+1;
end

%% GOOD: Vertical bar plot w total and RP/DP and hist bar
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.005],[.08 .08],[.08 .08]); %[vgap hgap], marg_h -[low up],marg_w -[l r]
%Matlab doesnt allow double axis plots to be rotated w the cam option so manually rotate elsewhere
% Add zero cell to create gaps
subPot_GWh(9,26)=0;
sel2plot=[1,26,2:13,26,14:25];

figure
colororder(c_ssp26(sel2plot,:))
% First plot theoretical
subplot(6,1,6)
% plot mid and far terms as separate
bar(diag(subPot_GWh(9,sel2plot)/1000'),'stacked','EdgeAlpha',0)
yline(subPot_GWh(9,1)/1000,'color','black','LineWidth',1.5,'LineStyle',':')
xline([2,15],'color',mygraylines,'LineWidth',1.05) % tf splits
set(gca,'YDir','reverse','YColor',[0.00,0.45,0.74],'XTick',1:25,'XTickLabels',[]) %,...
    %'XAxisLocation','top') %,'XLim',[0 110]),'YDir','reverse',
ylabel(pottypes5(1)) % "TWh per yr"],'FontWeight','bold')
ylim([0 2500])
text([8.5 15+6.5],[1 1]*2800,strcat(tname,": ", tframe),'HorizontalAlignment','center','FontWeight','bold')%,'Rotation',90)
applymyplotformat('')

% Then plot other pots
pp=4;
for pidx=[1,2,4]
    pt=newpots4(pidx);
    selp=find(tot_allscen.pottype==pt);
    % add gap
    selrows=[97; selp(1:12); 97; selp(13:24)];
    % plot total for hist and future
    subplot(6,1,pp)
    idata=[histHP.("Total TWh/yr")(histHP.pottype==pt);   tot_allscen.("Total TWh/yr")(selrows)]; % add historical
    b1=bar(diag(idata),'stacked','EdgeAlpha',0);
    hold all
    % add boxes for RP/DP totals for hist and future
    idata=[histHP{(histHP.pottype==pt), 2:3};   tot_allscen{selrows, 2:3}]; % add historical
    b2=bar(idata,.7,  'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.5,'LineStyle','-');
    ylabel(pottypes3(min(pidx,3))) %,'FontWeight','bold' )
    
    % set(gca,    'YDir','reverse')
    % Change color of  bar
    for k = 1:2
        b2(k).EdgeColor = mycolors{k}*0.6;
    end

    % add lines for hist pot and separaters for tfs
    b3=yline(tot_allscen.("Historical TWh/yr")(selp(1)),'color','black','LineWidth',1.5,'LineStyle',':');
    xticks(1:(length(selrows)+1))
    xline([2,15],'color',mygraylines,'LineWidth',1.05) % tf splits
    ylim([0 440])
    if pidx==2
        text(-1,210,"Energy in TWh per year",'HorizontalAlignment','center','FontWeight','bold','Rotation',90)
    end

    % add diamonds for RP/DP hist and then fut
    yyaxis right
    scshift=0;
    idata=[histHP{(histHP.pottype==pt), 5}; tot_allscen{selrows, 5}];
    s1(1)=scatter([1:length(selrows)+1]-scshift,idata,30,cl_RP*.6,'d','filled');
    idata=[histHP{(histHP.pottype==pt), 5}+histHP{(histHP.pottype==pt), 6}; tot_allscen{selrows, 6}+tot_allscen{selrows, 5}];
    s1(2)=scatter([1:length(selrows)+1]-scshift,idata,30, cl_DP*.6,'d','filled'); %manually stack scatter
    
    if pidx==2
        ylabel('Number of projects','FontWeight','bold')
    end
    ax = gca;
    ax.YAxis(2).Color = extraaxis ;
    %ylim([0 15000])

    applymyplotformat('')

    % rotate subplot
    % camroll(-90)

    % add legend and xticklabels
    if pp==4
        xticklabels([histrcpcornernames(1) '' rcpcornernames '' rcpcornernames])
        xtickangle(90)                    
    else
        xticklabels([])
    end

    pp=pp-1;
end
l=legend([b1(3:14) b2 s1 b1(1)],[repelem({'',},12-4) cornernames  'RP energy' 'DP energy' 'RP number' 'DP number' 'Historical'],'NumColumns',5, 'Location','northoutside');
l.Title.String=strjoin(rcpnames,'      ');
% rotate subplot so plot becomes horizontal bar
%camroll(-90)

%% Compile % change in different potential
% compile reshaped pots
compilepot_dHPpot=table([rcpcornernames'; rcpcornernames']);
compilepot_dHPpot.Theoretical=subPot_prctchange(9,2:25)';
compilepot_dHPnplants=table([rcpcornernames'; rcpcornernames']);

pp=1;
for pt=newpots4([1,2,4])
    selp=find(tot_allscen.pottype==pt);

    % reshaped pots by pot types for mid and far futures in one table
    compilepot_dHPpot.(pottypes3{pp})=changeprct_allscen.("Total TWh/yr")(selp);
    compilepot_dHPnplants.(pottypes3{pp})=changeprct_allscen.("Total number")(selp);
    pp=pp+1;
end


%% Compile total pot different potential
% compile reshaped pots
compilepot_TWhYrHPpot=table(histrcpcornernames','VariableNames',{'TWHpot'});
compilepot_TWhYrHPpot.Theoretical=subPot_GWh(9,1:25)'/1000;

pp=1;
for pt=newpots4([1,2,4])
    idata=[histHP.("Total TWh/yr")(histHP.pottype==pt);   tot_allscen.("Total TWh/yr")(tot_allscen.pottype==pt)]; % add historical
    % reshaped pots by pot types for mid and far futures in one table
    compilepot_TWhYrHPpot.(pottypes3{pp})=idata;
    pp=pp+1;
end


%% Bubble chart dP vs dT vs dHP
%dP_dT=readtable(fullfile(rootof,"FutRuns_Figs_Analysis" ,"FutScen_Tracker.xlsx"),Sheet="GCMdetails", Range="A1:H13");

figure;
for pp=2:5
    subplot(1,2,1)
    hold on
    bubblechart(dP_dT.x_P___,dP_dT.x_T_K_,compilepot_dHP{1:12,pp},'MarkerFaceAlpha',0);
    subplot(1,2,2)
    hold on
    bubblechart(dP_dT.x_P___,dP_dT.x_T_K_,compilepot_dHP{13:24,pp},'MarkerFaceAlpha',0);
end
for sp=1:2
    subplot(1,2,sp)
    ylabel("\Delta Change in temperature")
    xlabel("% Change in precipitation")
    bubblelim([-2 60])
    bubblesize([5 40])
    applymyplotformat(tf_full(sp))
end
bubblelegend("% Change in HP")
legend(pottypes5{1:4},'Orientation','Horizontal')

%% Write totals to matlab and excel
if writeSummary2file
    % Save tech, fin, sust pot in cleaned up and reordered version for
    % future and historical in one file
    save(matfile_summary,'tot_allscen','histHP','potlist')

    % Write the all scenarios compiled data
    tmp=sprintf("Last updated on: %s",datestr(datetime('now')));
    tot_allscen(size(tot_allscen,1)+1,1)={tmp{:}};
    writetable(tot_allscen,xlsfile_summary,'Sheet',"Compiled");
    writetable(histHP,xlsfile_summary,'Sheet',"Historical");

    writetable(change_allscen,xlsfile_summary,'Sheet',"Change_mag");
    writetable(changeprct_allscen,xlsfile_summary,'Sheet',"Change_prctofhisttot");
%rearranged
    writetable(compilepot_TWhYrHPpot,xlsfile_summary,'Sheet',"PotTWh_yr_rearranged");

    writetable(compilepot_dHPpot,xlsfile_summary,'Sheet',"PotTWh_yr_rearranged%change");
writetable(compilepot_dHPnplants,xlsfile_summary,'Sheet',"PotNplants_rearranged%change");

    disp("Written to file")
end