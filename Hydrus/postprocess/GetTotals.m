% Post-process to get total energy and # of projects for various scenario
% runs. Loops through folders based on folder name pattern specified
% option to write separate excel files for CAS
% Save basin level total Tech, Fin, Econ potential for all scenarios. For
% SAs addtionally saves delta change vals for Q and other to mat and only others to excel.

% Creates generic bar plots for all scenarios
% Creates cleaned up bar plots for 3 energy, 2 geohazard and multiple plot types for SA cases

clc
clearvars
close all
addpath('Hydrus\postprocess') % Add folder path
run('myVarNames.m')


plotmainfigs=0;
plotallSAfigs=0;
writeSummary2file=0;
ofname='Compile_Scenario_Totals03';
xlsfile_summary=fullfile(rootf,'output',[ofname '.xlsx']);
matfile_summary=fullfile(rootf,'output',[ofname '.mat']);
printsummary=0;


write4CAS=0;
xlspath=fullfile(rootf,'output','CAS_19July2022'); if ~isfolder(xlspath); mkdir(xlspath); end
xlsVarnames={'pid', 'Latss', 'Lonss', 'SysIDs', 'COEAlls_[USD2010/kWh]', 'PnetAlls_[GWh/yr]'};

idx=1;
addpath(genpath(fullfile(rootf,"Hydrus")))
selmixsust=9;

%% Loop through and compile FULL + REMAIN 3 policy types
runnameidx=103;
runname_prefix='R103_Energy'; %sprintf('R%d_Energy',runnameidx);
for scen={'Full', 'Remain'}%   'Full';  % 'Remain' ;  %
    for policyname={'Large','Medium','Mixed'}
        %% Separate tech and econ data
        scenario=scen{1};
        runname=strjoin({runname_prefix,scenario, policyname{1}},'_');

        % Prepare exact filename for matfile
        load(fullfile(rootf, 'output',scenario, continent_in,[runname '_Tech_Fin'],...
            sprintf('COEPOT_b%d_do.mat', nbasin)))

        %[tech_outdata_ss fin_outdata_ss] = getTechFinOutput(fname);
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
            fprintf('%s Financial Potential: %0.2f GWh\n',scenario, nansum(PnetAlls(COEAlls<=0.1)))
            fprintf('# of selected projects: %0.0f out of %0.0f\n',sum(~isnan(PnetAlls(COEAlls<=0.1))),length(PnetAlls))
            fprintf('# of selected River type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==2 & COEAlls<=0.1),sum(PnetAlls(SysIDs==2 & COEAlls<=0.1)))
            fprintf('# of selected Diversion type projects: %0.0f with %0.2f GWh\n\n',sum(SysIDs==1 & COEAlls<=0.1),sum(PnetAlls(SysIDs==1 & COEAlls<=0.1)))
        end

        % Evaluate % of projects below <=0.1
        prctfin(idx) = sum(~isnan(PnetAlls(COEAlls<=0.1)))/sum(~isnan(PnetAlls))*100;

        %Archive
        compile_PNetalls(idx,:)= [sum(PnetAlls(SysIDs==2)) sum(PnetAlls(SysIDs==1)) sum(SysIDs==2) sum(SysIDs==1)]; % RP/DP GWh/yr and counts
        compile_PNetalls(idx+1,:)= [sum(PnetAlls(SysIDs==2 & COEAlls<=0.1)) sum(PnetAlls(SysIDs==1 & COEAlls<=0.1)) sum(SysIDs==2 & COEAlls<=0.1) sum(SysIDs==1 & COEAlls<=0.1)]; %Fin - RP/DP GWh/yr and counts
        compile_runname{idx}=[runname '_Tech'];
        compile_runname{idx+1}= [runname '_Fin'];

        %% Get sust data
        load(fullfile(rootf, 'output',scenario, continent_in,[runname '_Sust_' geohazard_scennames{3}],sprintf('COEPOT_b%d_do.mat', nbasin)))
        selecon= COEAlls<=costlim;
        sust_outdata_ss=table(PID(selecon),latss(selecon),lonss(selecon),SysIDs(selecon),COEAlls(selecon),PnetAlls(selecon),...); % exported to excel
            'VariableNames',xlsVarnames);
        if printsummary
            fprintf('%s Sustainable Potential: %0.2f GWh\n',scenario, nansum(PnetAlls(COEAlls<=0.1)))
            fprintf('# of selected projects: %0.0f out of %0.0f\n',sum(~isnan(PnetAlls(COEAlls<=0.1))),length(PnetAlls))
            fprintf('# of selected River type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==2 & COEAlls<=0.1),sum(PnetAlls(SysIDs==2 & COEAlls<=0.1)))
            fprintf('# of selected Diversion type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==1 & COEAlls<=0.1),sum(PnetAlls(SysIDs==1 & COEAlls<=0.1)))
        end

        % Evaluate % of projects below <=0.1
        prctfin(idx+2) = sum(~isnan(PnetAlls(COEAlls<=0.1)))/sum(~isnan(PnetAlls))*100;

        %Archive
        compile_PNetalls(idx+2,:)= [sum(PnetAlls(SysIDs==2 & COEAlls<=0.1)) sum(PnetAlls(SysIDs==1 & COEAlls<=0.1)) sum(SysIDs==2 & COEAlls<=0.1) sum(SysIDs==1 & COEAlls<=0.1)]; % RP/DP GWh/yr and counts
        compile_runname{idx+2}= [runname '_Sust'];
        idx=idx+3;

        %% Write Tech/Fin/Sust to excel
        if write4CAS
            oxlsfile=fullfile(xlspath,[runname '.xlsx']);
            writetable(tech_outdata_ss,oxlsfile,'Sheet','Tech_RiskAverse')
            writetable(fin_outdata_ss,oxlsfile,'Sheet','Fin_RiskAverse')
            writetable(sust_outdata_ss,oxlsfile,'Sheet','Sust_RiskAverse')
        end
    end
end
%end
%% Add 3 GeoHazard Type, 2 GeoHaz rep, 3 Tech2Sust cases for Full_Mixed_Sust case
scenario='Full';
policyname='Mixed';
%runname_prefix={'R101_HazTyp', 'R101_HazRep'};
for runname_prefix={'R103_HazTyp', 'R103_HazRep', 'R103_Tech2Sust'}
    runname=strjoin({runname_prefix{1},scenario, policyname},'_');
    oxlsfile=fullfile(xlspath,[runname '.xlsx']);

    runfiles =  path2fldrfiles(fullfile(rootf, 'output',scenario, continent_in),...
        strcat(runname_prefix{1},"*"));

    tmp=split(runfiles,runname_prefix{1});
    geohazard_suffix=tmp(:,2);

    for g=1:length(runfiles)
        load(fullfile(runfiles{g},sprintf('COEPOT_b%d_do.mat', nbasin)))
        %PID=transpose(1:numel(latss));
        selecon= COEAlls<=costlim;
        sust_outdata_ss=table(PID(selecon),latss(selecon),lonss(selecon),SysIDs(selecon),COEAlls(selecon),PnetAlls(selecon),...); % exported to excel
            'VariableNames',xlsVarnames);

        if printsummary
            disp([runname geohazard_suffix{g}])
            fprintf('%s Sustainable Potential: %0.2f GWh\n',scenario, nansum(PnetAlls(COEAlls<=0.1)))
            fprintf('# of selected projects: %0.0f out of %0.0f\n',sum(~isnan(PnetAlls(COEAlls<=0.1))),length(PnetAlls))
            fprintf('# of selected River type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==2 & COEAlls<=0.1),sum(PnetAlls(SysIDs==2 & COEAlls<=0.1)))
            fprintf('# of selected Diversion type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==1 & COEAlls<=0.1),sum(PnetAlls(SysIDs==1 & COEAlls<=0.1)))
        end
        if write4CAS
            xx=geohazard_suffix{g};
            writetable(sust_outdata_ss,oxlsfile,'Sheet',xx(1:20));
        end

        compile_PNetalls(idx,:)= [sum(PnetAlls(SysIDs==2 & COEAlls<=0.1)) sum(PnetAlls(SysIDs==1 & COEAlls<=0.1)) sum(SysIDs==2 & COEAlls<=0.1) sum(SysIDs==1 & COEAlls<=0.1)]; % RP/DP GWh/yr and counts
        compile_runname{idx}= extractAfter(runfiles{g},[continent_in,'\']);
        idx=idx+1;
    end
end

mainidx_end=idx-1;

%% Add SA runs for Sust
scenario='Full';
policyname='Mixed';
ii=1;
for runname_prefix={'R103_SA_R', 'R103_SA_M'}  % first the SA_Q runs and then the SA_Other runs
    runname=strjoin({runname_prefix{1},scenario, policyname},'_');
    oxlsfile=fullfile(xlspath,[runname '.xlsx']);

    runfiles =  path2fldrfiles(fullfile(rootf, 'output',scenario, continent_in),...
        strcat(runname_prefix{1},"*Sust*"));
    for g=1:length(runfiles)
        load(fullfile(runfiles{g},sprintf('COEPOT_b%d_do.mat', nbasin)))
        %PID=transpose(1:numel(latss));
        selecon= COEAlls<=costlim;
        sust_outdata_ss=table(PID(selecon),latss(selecon),lonss(selecon),SysIDs(selecon),COEAlls(selecon),PnetAlls(selecon),...); % exported to excel
            'VariableNames',xlsVarnames);

        if printsummary
            disp([runname '_Sust_' geohazard_scennames{g}])
            fprintf('%s Sustainable Potential: %0.2f GWh\n',scenario, nansum(PnetAlls(COEAlls<=0.1)))
            fprintf('# of selected projects: %0.0f out of %0.0f\n',sum(~isnan(PnetAlls(COEAlls<=0.1))),length(PnetAlls))
            fprintf('# of selected River type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==2 & COEAlls<=0.1),sum(PnetAlls(SysIDs==2 & COEAlls<=0.1)))
            fprintf('# of selected Diversion type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==1 & COEAlls<=0.1),sum(PnetAlls(SysIDs==1 & COEAlls<=0.1)))
        end
        if write4CAS
            writetable(sust_outdata_ss,oxlsfile,'Sheet',['Sust_' geohazard_scennames{g}]);
        end

        compile_PNetalls(idx,:)= [sum(PnetAlls(SysIDs==2 & COEAlls<=0.1)) sum(PnetAlls(SysIDs==1 & COEAlls<=0.1)) sum(SysIDs==2 & COEAlls<=0.1) sum(SysIDs==1 & COEAlls<=0.1)]; % RP/DP GWh/yr and counts
        compile_runname{idx}= extractAfter(runfiles{g},[continent_in,'\']);
        idx=idx+1;
    end
    SA_Sust_idx_end(ii)=idx-1;ii=ii+1;
end

%% Add SA runs for Tech
scenario='Full';
policyname='Mixed';
ii=1;
for runname_prefix={'R103_SA_M'}
    runname=strjoin({runname_prefix{:},scenario, policyname},'_');
    oxlsfile=fullfile(xlspath,[runname '.xlsx']);

    runfiles =  path2fldrfiles(fullfile(rootf, 'output',scenario, continent_in),...
        strcat(runname_prefix{:},"*Tech_Fin*"));
    for g=1:length(runfiles)
        load(fullfile(runfiles{g},sprintf('COEPOT_b%d_do.mat', nbasin)))
        %PID=transpose(1:numel(latss));
        selall= ~isnan(COEAlls);
        sust_outdata_ss=table(PID(selall),latss(selall),lonss(selall),SysIDs(selall),COEAlls(selall),PnetAlls(selall),...); % exported to excel
            'VariableNames',xlsVarnames);

        if printsummary
            disp([runname '_Sust_' geohazard_scennames{g}])
            fprintf('%s Sustainable Potential: %0.2f GWh\n',scenario, nansum(PnetAlls(COEAlls<=0.1)))
            fprintf('# of selected projects: %0.0f out of %0.0f\n',sum(~isnan(PnetAlls(COEAlls<=0.1))),length(PnetAlls))
            fprintf('# of selected River type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==2 & COEAlls<=0.1),sum(PnetAlls(SysIDs==2 & COEAlls<=0.1)))
            fprintf('# of selected Diversion type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==1 & COEAlls<=0.1),sum(PnetAlls(SysIDs==1 & COEAlls<=0.1)))
        end
        if write4CAS
            writetable(sust_outdata_ss,oxlsfile,'Sheet',['Sust_' geohazard_scennames{g}]);
        end

        compile_PNetalls(idx,:)= [sum(PnetAlls(SysIDs==2)) sum(PnetAlls(SysIDs==1)) sum(SysIDs==2) sum(SysIDs==1)]; % RP/DP GWh/yr and counts
        compile_runname{idx}= strcat(extractAfter(runfiles{g},[continent_in,'\']),'_TECHPOT');
        idx=idx+1;
    end
end
SA_Tech_idx_end=idx-1;ii=ii+1;

%% Add SA runs for Fin
scenario='Full';
policyname='Mixed';
ii=1;
for runname_prefix={'R103_SA_M'}
    runname=strjoin({runname_prefix{1},scenario, policyname},'_');
    oxlsfile=fullfile(xlspath,[runname '.xlsx']);

    runfiles =  path2fldrfiles(fullfile(rootf, 'output',scenario, continent_in),...
        strcat(runname_prefix{1},"*Tech_Fin*"));
    for g=1:length(runfiles)
        load(fullfile(runfiles{g},sprintf('COEPOT_b%d_do.mat', nbasin)))
        %PID=transpose(1:numel(latss));
        selecon= COEAlls<=costlim;
        sust_outdata_ss=table(PID(selecon),latss(selecon),lonss(selecon),SysIDs(selecon),COEAlls(selecon),PnetAlls(selecon),...); % exported to excel
            'VariableNames',xlsVarnames);

        if printsummary
            disp([runname '_Sust_' geohazard_scennames{g}])
            fprintf('%s Sustainable Potential: %0.2f GWh\n',scenario, nansum(PnetAlls(COEAlls<=0.1)))
            fprintf('# of selected projects: %0.0f out of %0.0f\n',sum(~isnan(PnetAlls(COEAlls<=0.1))),length(PnetAlls))
            fprintf('# of selected River type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==2 & COEAlls<=0.1),sum(PnetAlls(SysIDs==2 & COEAlls<=0.1)))
            fprintf('# of selected Diversion type projects: %0.0f with %0.2f GWh\n',sum(SysIDs==1 & COEAlls<=0.1),sum(PnetAlls(SysIDs==1 & COEAlls<=0.1)))
        end
        if write4CAS
            writetable(sust_outdata_ss,oxlsfile,'Sheet',['Sust_' geohazard_scennames{g}]);
        end

        compile_PNetalls(idx,:)= [sum(PnetAlls(SysIDs==2 & COEAlls<=0.1)) sum(PnetAlls(SysIDs==1 & COEAlls<=0.1)) sum(SysIDs==2 & COEAlls<=0.1) sum(SysIDs==1 & COEAlls<=0.1)]; % RP/DP GWh/yr and counts
        compile_runname{idx}= strcat(extractAfter(runfiles{g},[continent_in,'\']),'_FINPOT');
        idx=idx+1;
    end
end
SA_Fin_idx_end=idx-1;ii=ii+1;

%% Compile all scenarios
cl=strrep(strrep(compile_runname,'R103_',''),'_','-')'; %Replace _ w - so names are not weird
%cl=categorical(cl,cl);
xlsVarnames={'pid', 'Latss', 'Lonss', 'SysIDs', 'COEAlls_[USD2010/kWh]', 'PnetAlls_[GWh/yr]'};
tot_allscen=table(cl, compile_PNetalls(:,1)/1000, compile_PNetalls(:,2)/1000, sum(compile_PNetalls(:,1:2),2)/1000,...
    compile_PNetalls(:,3), compile_PNetalls(:,4), sum(compile_PNetalls(:,3:4),2),...
    'VariableNames',{'Runname','TWh/yr of River power','TWh/yr of Diversion canal','Total TWh/yr','Number of River power','Number of Diversion canal','Total number'});
% reorder the SA_other vars
%xx=tot_allscen([1:53 54,57,56,58,60,59,55,61    64    63    65    67    66    62 68    71    70    72    74    73    69 75    78    77    79    81    80    76 82    85    84    86    88    87    83   89    92    91    93    95    94    90],:);
disp("Compiled all data")

all_prctfin_num=table(cl, prctfin','VariableNames',{'Runname','Prct of plants costing <=0.10'});

%% SA_Sust other: Get delta change for rotated bar and scatter
% SA_others_stacked_delta=zeros(14*3,2,2); % param ordered by tech/fin/sust x RP/DP x min/max
varnames={};
m=1;
pp=0;
for scen=1:14 %loop through SA param scens
    varnames=[varnames; cl([SA_Sust_idx_end(2)+scen SA_Tech_idx_end+scen SA_Sust_idx_end(1)+scen])];

    for plant=1:2 % RP/DP in separate columns
        % Eval delta for energy= tech/fin/sust_SA - tech/fin/sust_default %
        % % tot_allscen has Sust before Tech and Fin scenarios for SA other params
        SA_op_stacked_delta(pp+1,plant)=tot_allscen{SA_Sust_idx_end(2)+scen, 1+plant} - tot_allscen{7, 1+plant};%Tech
        SA_op_stacked_delta(pp+2,plant)=tot_allscen{SA_Tech_idx_end+scen, 1+plant} - tot_allscen{8, 1+plant};   %Fin
        SA_op_stacked_delta(pp+3,plant)=tot_allscen{SA_Sust_idx_end(1)+scen, 1+plant} - tot_allscen{9, 1+plant};%Sust

        SA_op_stacked_deltanum(pp+1,plant)=tot_allscen{SA_Sust_idx_end(2)+scen, 4+plant} - tot_allscen{7, 4+plant};%Tech
        SA_op_stacked_deltanum(pp+2,plant)=tot_allscen{SA_Tech_idx_end+scen, 4+plant} - tot_allscen{8, 4+plant};   %Fin
        SA_op_stacked_deltanum(pp+3,plant)=tot_allscen{SA_Sust_idx_end(1)+scen, 4+plant} - tot_allscen{9, 4+plant};%Sust

        % Eval delta as % of total = delta/(total RP+DP)*100
        SA_op_stacked_deltaprct(pp+1,plant)=SA_op_stacked_delta(pp+1,plant)/tot_allscen{7,4}*100;%Tech
        SA_op_stacked_deltaprct(pp+2,plant)=SA_op_stacked_delta(pp+2,plant)/tot_allscen{8, 4}*100;   %Fin
        SA_op_stacked_deltaprct(pp+3,plant)=SA_op_stacked_delta(pp+3,plant)/tot_allscen{9, 4}*100;%Sust


        SA_op_stacked_deltaprctnum(pp+1,plant)=SA_op_stacked_deltanum(pp+1,plant)/tot_allscen{7, 7}*100;%Tech
        SA_op_stacked_deltaprctnum(pp+2,plant)=SA_op_stacked_deltanum(pp+2,plant)/tot_allscen{8, 7}*100;   %Fin
        SA_op_stacked_deltaprctnum(pp+3,plant)=SA_op_stacked_deltanum(pp+3,plant)/tot_allscen{9, 7}*100;%Sust
    end
    pp=pp+3;
end

cl_varnames=erase(varnames,{'SA-','-Full-Mixed','-Tech-Fin','-RiskAverse'});
SA_Others_delta=table(cl_varnames,varnames, SA_op_stacked_delta, SA_op_stacked_deltanum, SA_op_stacked_deltaprct, SA_op_stacked_deltaprctnum,...
    'VariableNames',{'shortRunname','FullRunname','Delta TWh/yr RP/DP','Delta # RP/DP','Delta as % total TWh/yr RP/DP','Delta as % total # RP/DP'});


%% Main scenario figs
if plotmainfigs
    %% Plot bar for all main scenarios
    figure
    selidx=1:mainidx_end;
    b1=bar(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 2:3},'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',baralpha);
    grid minor
    ylabel('Bar for TWh per year')
    yyaxis right
    scatter(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 5},30,cl_RP*.6,'d','filled')
    hold on
    scatter(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 6}+tot_allscen{selidx, 5},30, cl_DP*.6,'d','filled') %manually stack scatter

    %Setup
    legend(tot_allscen.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
    ylabel('\diamondsuit for number of projects\diamondsuit ')
    xlabel('Scenarios')
    % Change color of  bar
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end

    %alpha(b1,0.7)
    ax = gca;
    ax.YAxis(2).Color = 'k';
    title("Main Scenarios")
    %% FINAL: Plot bar for 3 energy scenarios full and remain overlaid
    figure
    selfull=1:9; %mainidx_end;
    selremain=10:18; %mainidx_end;

    b1=bar(tot_allscen{selfull, 2:3},'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',0.7,'EdgeAlpha',0);
    hold all
    b2=bar(tot_allscen{selremain, 2:3},0.6, 'stacked','FaceColor',"flat",'EdgeColor',"flat",'FaceAlpha',0,'LineWidth',1.15,'LineStyle','-');
    ylim([0 350])

    % add lines separating boxes
    xline([3.5,6.5],'Color', 'k')


    % Change color of  bar
    ylabel('Bar = Energy in TWh/year','FontWeight','bold')

    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};

        b2(k).EdgeColor = mycolors{k}*0.6;
    end
    yyaxis right
    scshift=0.1;
    scatter([1:9]-scshift,tot_allscen{selfull, 5},30,cl_RP*.6,'d','filled')
    hold on
    scatter([1:9]-scshift,tot_allscen{selfull, 6}+tot_allscen{selfull, 5},30, cl_DP*.6,'d','filled') %manually stack scatter
    %
    scatter([1:9]+scshift,tot_allscen{selremain, 5},30,cl_RP*.6,'d','LineWidth',1.15)
    scatter([1:9]+scshift,tot_allscen{selremain, 6}+tot_allscen{selremain, 5},30, cl_DP*.6,'d','LineWidth',1.15) %manually stack scatter
    % fix yyaxis
    ylabel('\diamondsuit = Number of projects','FontWeight','bold')
    ax = gca;
    ax.YAxis(2).Color = 'k';

    % add scenario labels
    xticklabels(repmat([pottypes3_short],1,3))
    xlabel([""; "Energy Focus Scenarios"],'FontWeight','bold')
    text([2 5 8],-350*[1 1 1],searchtypes,'HorizontalAlignment','center','FontWeight','bold')
    % add legend
    l=legend(["Full", "Full", "Remain |", "Remain |","Full","Full", "Remain", "Remain"], "NumColumns", 4,'Location','northoutside'); %,'Orientation','horizontal')%)
    l.Title.String=sprintf("Energy Potential | Number of Projects ");

    % legend(["RP: Full Energy", "DP: Full Energy", "RP: Remaining Energy", "DP: Remaining Energy",...
    %     "RP: Full Number","DP: Full Number", "RP: Remaining Number","DP: Remaining Number"], "NumColumns", 4,'Location','northoutside') %,'Orientation','horizontal')%)
    grid minor
    % Create textbox
    annotation('textbox',...
        [0.232187761944677 0.855755894590846 0.575859178541493 0.0776699029126213], ...    %[0.184516444185945 0.852821668148887 0.215419501133787 0.0999999999999999],...
        'VerticalAlignment','bottom',...
        'String',['River power plant:',newline(),'Diversion canal plant:'],...
        'HorizontalAlignment','right',...
        'FontSize',14,...
        'FontName','segoe ui',...
        'FitBoxToText','on',...
        'EdgeColor',[0.149019607843137 0.149019607843137 0.149019607843137],...
        'BackgroundColor',[1 1 1]);


    %% FINAL: Plot bar for geo hazard type and geo hazard rep scenarios
    cl=tot_allscen.Runname;
    figure
    subplot(1,2,1)
    selidx=[19:21,9];
    b1=bar(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 2:3},'stacked','FaceColor',"flat",'EdgeColor',"flat");
    grid minor
    ylabel('Bar = Energy in TWh/year','FontWeight','bold')
    ylim([0 350])

    yyaxis right
    scatter(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 5},30,cl_RP*.6,'d','filled')
    hold on
    scatter(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 6}+tot_allscen{selidx, 5},30, cl_DP*.6,'d','filled') %manually stack scatter
    %Setup
    %legend(t.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
    %ylabel('\diamondsuit = Number of projects','FontWeight','bold')
    xticklabels([geohazard_names "All 3"])
    xlabel("Geo-hazard types",'FontWeight','bold')
    % Change color of  bar
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end

    alpha(b1,0.7)
    ax = gca;
    ax.YAxis(2).Color = 'k';


    subplot(1,2,2)
    selidx=[9,22:23];
    b1=bar(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 2:3},'stacked','FaceColor',"flat",'EdgeColor',"flat");
    grid minor
    ylim([0 350])

    %ylabel('Bar = Energy in TWh/year','FontWeight','bold')
    yyaxis right
    scatter(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 5},30,cl_RP*.6,'d','filled')
    hold on
    scatter(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 6}+tot_allscen{selidx, 5},30, cl_DP*.6,'d','filled') %manually stack scatter

    %Setup
    %legend(t.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
    ylabel('\diamondsuit = Number of projects','FontWeight','bold')
    xticklabels(geohazard_scennames_cl([3,2,4]))

    xlabel(['Geo-hazard risk', newline(), 'representations'],'FontWeight','bold')
    % Change color of  bar
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end

    alpha(b1,0.7)
    ax = gca;
    ax.YAxis(2).Color = 'k';

end
%% SA figs
if plotallSAfigs
    %% Plot bar for SA Q scenarios
    selidx=[9 mainidx_end+1:SA_Sust_idx_end(1)];
    figure; %subplot(3,1,2)
    varnames=extractBefore(cl(selidx),'-Full');
    varnames(1)={'Default'};

    b1=bar(categorical(varnames,varnames),tot_allscen{selidx, 2:3},'stacked');
    grid minor
    yyaxis right
    scatter(categorical(varnames,varnames),tot_allscen{selidx, 5},50,cl_RP*.6,'d','filled')
    hold on
    scatter(categorical(varnames,varnames),tot_allscen{selidx, 6}+tot_allscen{selidx, 5},50,cl_DP*.6,'d','filled') %manually stack scatter
    legend(tot_allscen.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end
    alpha(b1,0.7)
    ax = gca;
    ax.YAxis(2).Color = 'k';
    title("SA Qs Sust")

    %% Plot bar for SA_Sust other scenarios
    selidx=[9 SA_Sust_idx_end(1)+1:SA_Sust_idx_end(2)];
    figure; %subplot(3,1,3)
    varnames=extractBefore(cl(selidx),'-Full');
    varnames(1)={'Default'};

    b1=bar(categorical(varnames,varnames),tot_allscen{selidx, 2:3},'stacked');
    grid minor
    yyaxis right
    scatter(categorical(varnames,varnames),tot_allscen{selidx, 5},50,cl_RP*.6,'d','filled')
    hold on
    scatter(categorical(varnames,varnames),tot_allscen{selidx, 6}+tot_allscen{selidx, 5},50,cl_DP*.6,'d','filled') %manually stack scatter
    legend(tot_allscen.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end
    alpha(b1,0.7)
    ax = gca;
    ax.YAxis(2).Color = 'k';
    title("SA - OtherParms SUST")

    %% Plot bar for SA_Tech other scenarios
    selidx=[7 SA_Sust_idx_end(2)+1:SA_Tech_idx_end];
    figure; %subplot(3,1,3)
    varnames=extractBefore(cl(selidx),'-Full');%   cl(selidx); %
    varnames(1)={'Default - Tech'};

    b1=bar(categorical(varnames,varnames),tot_allscen{selidx, 2:3},'stacked');
    grid minor
    yyaxis right
    scatter(categorical(varnames,varnames),tot_allscen{selidx, 5},50,cl_RP*.6,'d','filled')
    hold on
    scatter(categorical(varnames,varnames),tot_allscen{selidx, 6}+tot_allscen{selidx, 5},50,cl_DP*.6,'d','filled') %manually stack scatter
    legend(tot_allscen.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end
    alpha(b1,0.7)
    ax = gca;
    ax.YAxis(2).Color = 'k';
    title("SA - OtherParms TECH")

    %% Plot bar for SA_Fin other scenarios
    selidx=[8 SA_Tech_idx_end+1:SA_Fin_idx_end];
    figure; %subplot(3,1,3)
    varnames=extractBefore(cl(selidx),'-Full');%   cl(selidx); %
    varnames(1)={'Default - FIN'};

    b1=bar(categorical(varnames,varnames),tot_allscen{selidx, 2:3},'stacked');
    grid minor
    yyaxis right
    scatter(categorical(varnames,varnames),tot_allscen{selidx, 5},50,cl_RP*.6,'d','filled')
    hold on
    scatter(categorical(varnames,varnames),tot_allscen{selidx, 6}+tot_allscen{selidx, 5},50,cl_DP*.6,'d','filled') %manually stack scatter
    legend(tot_allscen.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end
    alpha(b1,0.7)
    ax = gca;
    ax.YAxis(2).Color = 'k';
    title("SA - OtherParms FIN")

    %% FINAL SA_SustQ: Plot circle matrix of total val circles SA_SustQ
    RP_Qs =      {'Q30', 'Q40', 'Q50'};
    largeDP_Qs = {'Q25', 'Q30', 'Q40'};
    smallDP_Qs = {'Q70', 'Q80', 'Q90'};

    % Rescale vals for scatter plot sizing
    selidx=mainidx_end+1:SA_Sust_idx_end(1);
    scattervals=rescale(tot_allscen{selidx, 2:3},20,2000);
    runcnt=1;

    figure;
    for ii=1:numel(RP_Qs)
        subplot(1,3,ii); title(sprintf('RP = %s',RP_Qs{ii}))
        hold all
        for jj=1:numel(largeDP_Qs)
            for kk=1:numel(smallDP_Qs)
                scatter(kk,jj,scattervals(runcnt,2)+scattervals(runcnt,1),cl_RP,'filled','MarkerFaceAlpha',0.7) %RP
                scatter(kk,jj,scattervals(runcnt,2),cl_DP,'filled','MarkerFaceAlpha',0.7) %DP
                runcnt=runcnt+1;
            end
        end

        grid on
        axis square
        box on
        set(gca,'Xtick',[1:3],'XTickLabel',smallDP_Qs,'Ytick',[1:3],'YTickLabel',largeDP_Qs)
        xlim([0 4]);ylim([0 4])
        xlabel('smallDP_Qs')
        ylabel('largeDP_Qs')
    end
    legend(planttypes)
    sgtitle("SA Qs total energy magnitude")

    %% SA_SustQ: Get delta change
    selidx=mainidx_end+1:SA_Sust_idx_end(1);
    SA_Q_delta=tot_allscen{selidx, 2:3}-tot_allscen{selmixsust,2:3};
    SA_Q_delta_prct=SA_Q_delta./tot_allscen{selmixsust,2:3}*100;

    SA_Q_delta_num=tot_allscen{selidx, 5:6}-tot_allscen{selmixsust,5:6};
    SA_Q_delta_num_prct=SA_Q_delta_num./tot_allscen{selmixsust,5:6}*100;
    SA_Q_varnames=extractBefore(cl(selidx),'-Full');

    %% SA_SustQ: BAD Plot matrix of delta val circles SA_SustQ - DOES NOT SHOW NEGATIVES
    RP_Qs =      {'Q30', 'Q40', 'Q50'};
    largeDP_Qs = {'Q25', 'Q30', 'Q40'};
    smallDP_Qs = {'Q70', 'Q80', 'Q90'};

    % Rescale vals for scatter plot sizing
    scattervals_delta=rescale(SA_Q_delta,20,2000);
    runcnt=1;

    figure;
    for ii=1:numel(RP_Qs)
        subplot(1,3,ii); title(sprintf('RP = %s',RP_Qs{ii}))
        hold all
        for jj=1:numel(largeDP_Qs)
            for kk=1:numel(smallDP_Qs)
                scatter(kk,jj,scattervals_delta(runcnt,2)+scattervals_delta(runcnt,1),cl_RP,'filled','MarkerFaceAlpha',0.7) %RP
                scatter(kk,jj,scattervals_delta(runcnt,2),cl_DP,'filled','MarkerFaceAlpha',0.7) %DP
                runcnt=runcnt+1;
            end
        end

        grid on
        axis square
        box on
        set(gca,'Xtick',[1:3],'XTickLabel',smallDP_Qs,'Ytick',[1:3],'YTickLabel',largeDP_Qs)
        xlim([0 4]);ylim([0 4])
        xlabel('smallDP_Qs')
        ylabel('largeDP_Qs')
    end
    legend(planttypes)
    sgtitle("SA Qs change in energy ")

    %% SA_SustQ: Plot simple bar for delta SA Q scenarios
    selidx=[mainidx_end+1:SA_Sust_idx_end(1)];
    figure; %subplot(3,1,2)
    varnames=extractBefore(cl(selidx),'-Full');

    b1=bar(categorical(varnames,varnames),tot_allscen{selidx, 2:3}-tot_allscen{selmixsust,2:3},'stacked');
    grid minor
    yyaxis right
    scatter(categorical(varnames,varnames),tot_allscen{selidx, 5}-tot_allscen{selmixsust,5},50,cl_RP*.6,'d','filled')
    hold on
    scatter(categorical(varnames,varnames),tot_allscen{selidx, 6}-tot_allscen{selmixsust,6}+tot_allscen{selidx, 5}-tot_allscen{selmixsust,5},50,cl_DP*.6,'d','filled') %manually stack scatter
    legend(tot_allscen.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end
    alpha(b1,0.7)
    ax = gca;
    ax.YAxis(2).Color = 'k';
    title("Delta change for SA Qs Sust")

    %% FINAL SA_SustQ: Plot bar matrix for delta SA Q scenarios
    tt=1;
    varnames=extractBefore(cl,'-Full');
    runcnt=mainidx_end+1;
    f=figure; %subplot(3,1,2)
    for ii=1:numel(RP_Qs)
        hold all
        for jj=1:numel(largeDP_Qs)
            subtightplot (3, 3, tt, [0.025 0.08],[.05 .1],[.1 .06]); %[vgap hgap], marg_h,marg_w subplot(3,3,tt);
            selidx=runcnt:runcnt+2;
            b1=bar(categorical(varnames(selidx),varnames(selidx)),tot_allscen{selidx, 2:3}-tot_allscen{selmixsust,2:3},'stacked');
            for k = 1:2
                b1(k).FaceColor = mycolors{k};
                b1(k).EdgeColor = mycolors{k};
            end
            alpha(b1,0.7)
            ylim([-22 22])
            if jj==1 % add left ylabel to left col
                ylabel('\Delta TWh/yr')
                text(-0.8,-5.5,sprintf('RP = %s',RP_Qs{ii}),'FontWeight','bold','Rotation',90,'FontSize',12);%
            end

            yyaxis right
            scatter(categorical(varnames(selidx),varnames(selidx)),(tot_allscen{selidx, 5}-tot_allscen{selmixsust,5}),50,cl_RP*.6,'d','filled')
            hold on
            scatter(categorical(varnames(selidx),varnames(selidx)),(tot_allscen{selidx, 6}-tot_allscen{selmixsust,6})+(tot_allscen{selidx, 5}-tot_allscen{selmixsust,5}),50,cl_DP*.6,'d','filled') %manually stack scatter
            ylim([-110 110])
            if jj==3 % add right ylabel to right col
                ylabel('\Delta # of plants')
            end

            ax = gca;
            ax.YAxis(2).Color = 'k';
            set(gca,'XTickLabel',[])
            grid minor
            box on
            if ii==2 && jj==2
                % Add rectangle to show default
                text(1.65,100,"Default",'FontAngle','italic')
                rectangle('Position',[1.5 -108 1 217],'LineStyle','-.')
            end
            if ii==1  % add title to top row
                title(sprintf('Large DP = %s',largeDP_Qs{jj}),'FontSize',12);
            elseif ii==3  % add xlabel to bottom row
                set(gca,'XTickLabel',smallDP_Qs)
                text(1.6, -135,'Small DP','FontWeight','bold','FontSize',12);%
            end
            tt=tt+1;
            runcnt=runcnt+3;
        end
    end
    sgtitle("SA Qs Sust - Delta change")
    legend(tot_allscen.Properties.VariableNames{[2:3,5:6]}) %,'Location','southoutside','Orientation','horizontal')%)

    %% FINAL SA_SustQ: Plot bar matrix for % delta SA Q scenarios
    tt=1;
    runcnt=1;
    f=figure; %subplot(3,1,2)
    for ii=1:numel(RP_Qs)
        hold all
        for jj=1:numel(largeDP_Qs)
            subtightplot (3, 3, tt, [0.025 0.08],[.08 .11],[.1 .1]); %[vgap hgap], marg_h,marg_w subplot(3,3,tt);
            selidx=runcnt:runcnt+2;
            b1=bar(categorical(SA_Q_varnames(selidx),SA_Q_varnames(selidx)),SA_Q_delta_prct(selidx,:),'stacked');
            for k = 1:2
                b1(k).FaceColor = mycolors{k};
                b1(k).EdgeColor = mycolors{k};
            end
            alpha(b1,0.7)
            ylim([-45 35])
            if jj==1 % add left ylabel to left col
                if ii==2; ylabel('% change in TWh/yr','FontSize',12); end
                text(-1.5,-5.5,sprintf('RP = %s',RP_Qs{ii}),'FontWeight','bold','Rotation',90,'FontSize',12,'HorizontalAlignment','center');%
            end

            yyaxis right
            scatter(categorical(SA_Q_varnames(selidx),SA_Q_varnames(selidx)),SA_Q_delta_num_prct(selidx,1),50,cl_RP*.6,'d','filled','MarkerFaceAlpha',baralpha)
            hold on
            scatter(categorical(SA_Q_varnames(selidx),SA_Q_varnames(selidx)),SA_Q_delta_num_prct(selidx,1)+SA_Q_delta_num_prct(selidx,2),50,cl_DP*.6,'d','filled','MarkerFaceAlpha',baralpha) %manually stack scatter
            ylim([-45 35])
            if jj==3 % add right ylabel to right col
                if ii==2;  ylabel('% change in # of plants','FontSize',12); end
            end

            ax = gca;
            ax.YAxis(2).Color = 'k';
            set(gca,'XTickLabel',[])
            grid minor
            box on
            if ii==2 && jj==2
                % Add rectangle to show default
                text(2,31,"Default",'FontAngle','italic','HorizontalAlignment','center')
                rectangle('Position',[1.5 -44 1 34+44],'LineStyle','-.')
            end
            if ii==1  % add title to top row
                title(sprintf('Large DP = %s',largeDP_Qs{jj}),'FontSize',12);
            elseif ii==3  % add xlabel to bottom row
                set(gca,'XTickLabel',smallDP_Qs)
                if jj==2; xlabel('Small DP','FontWeight','bold','FontSize',12); end
            end
            tt=tt+1;
            runcnt=runcnt+3;
        end
    end
    sgtitle("SA Qs Sust - Delta change in %")
    legend(tot_allscen.Properties.VariableNames{[2:3,5:6]},'NumColumns',2) %,'Location','southoutside','Orientation','horizontal')%)


    %% SA_Sust other: Plot Rotated bar for SA_Sust other params
    figure
    selidx= [9 SA_Sust_idx_end(1)+1+7:SA_Sust_idx_end(2)]; %Min
    subplot(1,2,1)
    b1=barh(categorical(cl(selidx),cl(selidx)),tot_allscen{selidx, 2:3},'stacked')
    set(gca,'XDir','reverse','YAxisLocation','right','YDir','reverse','YTickLabels',[],'XLim',[0 110])
    grid on
    grid minor
    title("Lower")

    selidx= [9 SA_Sust_idx_end(1)+1:SA_Sust_idx_end(1)+7];
    varnamesCL=[".            DEFAULT            .";
        ".          minH-smallDP         .";
        ".generation efficiency - largeDP.";
        ".generation efficiency - smallDP.";
        ".generation efficiency - largeRD.";
        ".         interest rate         .";
        ".      environmental flow       .";
        ".    reservoir storage limit    .";];  %extractBetween(cl(SA_idx_end(1)+1:SA_idx_end(1)+7),"Max-","-Full");

    subplot(1,2,2)
    b2=barh(categorical(varnamesCL,varnamesCL),tot_allscen{selidx, 2:3},'stacked');%,'FaceAlpha',0.5)
    set(gca,'YDir','reverse')%,'YTickLabels',varnamesCL)
    xlim([0 110])
    grid on
    grid minor
    title("Higher")
    legend(tot_allscen.Properties.VariableNames(2:3),'Location','bestoutside')


    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
        b2(k).FaceColor = mycolors{k};
        b2(k).EdgeColor = mycolors{k};
    end
    sgtitle("SA Sust OtherParms")

    %b(1).CData(1,:) = [.5 0 .5]; %only first bar



    %% SA_Sust other: bad Grouped stacked bars - absolute
    SA_others_stacked=zeros(14,3,2); % param x tech/fin/sust x RP/DP
    for p=1:14
        for plant=1:2
            SA_others_stacked(p,1,plant)=tot_allscen{SA_Sust_idx_end(2)+p, 1+plant};
            SA_others_stacked(p,2,plant)=tot_allscen{SA_Tech_idx_end+p, 1+plant};
            SA_others_stacked(p,3,plant)=tot_allscen{SA_Sust_idx_end(1)+p, 1+plant};
        end
    end
    planttypes={'River Power','Diversion Canal'};
    stackLabels={'Tech', 'Fin', 'Sust'};
    groupLabels= extractBefore(cl(SA_Sust_idx_end(1)+1:SA_Sust_idx_end(1)+14),'-Full');%   cl(selidx); %

    plotBarStackGroupsDualYaxis(SA_others_stacked, ...
        'groupLabels', groupLabels, ...
        'stackLabels',stackLabels, ...
        'legendEntrys', planttypes)

    % selidx=[8 SA_Fin_idx_end(1)+1:numel(cl)];
    % selidx=[7 SA_Tech_idx_end(1)+1:SA_Tech_idx_end(2)];
    % selidx=[9 SA_Sust_idx_end(1)+1:SA_Sust_idx_end(2)];

    %% SA_Sust other: bad Get delta change
    SA_others_stacked_delta=zeros(14,3,2); % param x tech/fin/sust x RP/DP
    for p=1:14
        for plant=1:2
            % Eval delta= tech/fin/sust_SA - tech/fin/sust_default
            SA_others_stacked_delta(p,1,plant)=tot_allscen{SA_Sust_idx_end(2)+p, 1+plant} - tot_allscen{7, 1+plant};%Tech
            SA_others_stacked_delta(p,2,plant)=tot_allscen{SA_Tech_idx_end+p, 1+plant} - tot_allscen{8, 1+plant};   %Fin
            SA_others_stacked_delta(p,3,plant)=tot_allscen{SA_Sust_idx_end(1)+p, 1+plant} - tot_allscen{9, 1+plant};%Sust

            SA_others_stacked_deltaprct(p,1,plant)=SA_others_stacked_delta(p,1,plant)/tot_allscen{7, 1+plant};%Tech
            SA_others_stacked_deltaprct(p,2,plant)=SA_others_stacked_delta(p,2,plant)/tot_allscen{8, 1+plant};   %Fin
            SA_others_stacked_deltaprct(p,3,plant)=SA_others_stacked_delta(p,3,plant)/tot_allscen{9, 1+plant};%Sust
        end
    end
    % SA_Sust other: BAD grouped stacked bars This ignores negative vals so not useful really
    %figure
    plotBarStackGroupsDualYaxis(SA_others_stacked_delta(1:7,:,:), ...
        'groupLabels', groupLabels(1:7), ...
        'stackLabels',stackLabels, ...
        'legendEntrys', planttypes)
    title('Increase')
    plotBarStackGroupsDualYaxis(SA_others_stacked_delta(8:14,:,:), ...
        'groupLabels', groupLabels(8:14), ...
        'stackLabels',stackLabels, ...
        'legendEntrys', planttypes)
    title('Decrease')

    %% SA_Sust other: Simple bar plot delta to see -ves
    varnames=groupLabels; %
    figure;
    subplot(2,1,1)
    b1=bar(categorical(varnames,varnames),SA_others_stacked_delta(:,:,1));
    title("RP")
    grid minor
    legend(stackLabels)%,'Location','northoutside','Orientation','horizontal')%)
    alpha(b1,0.7)

    subplot(2,1,2)
    b1=bar(categorical(varnames,varnames),SA_others_stacked_delta(:,:,2));
    title("DP")
    grid minor
    legend(stackLabels)%,'Location','northoutside','Orientation','horizontal')%)
    alpha(b1,0.7)
    sgtitle("SA - diff for other params")

    %% SA_Sust other: Simple bar plot delta PRCT to see -ves
    varnames=groupLabels; %
    figure;
    subplot(2,1,1)
    b1=bar(categorical(varnames,varnames),SA_others_stacked_deltaprct(:,:,1));
    title("RP")
    grid minor
    legend(stackLabels)%,'Location','northoutside','Orientation','horizontal')%)
    alpha(b1,0.7)

    subplot(2,1,2)
    b1=bar(categorical(varnames,varnames),SA_others_stacked_deltaprct(:,:,2));
    title("DP")
    grid minor
    legend(stackLabels)%,'Location','northoutside','Orientation','horizontal')%)
    alpha(b1,0.7)
    sgtitle("SA - diff prct for other params")

    %% SA_Sust other: basic bar plot of delta change in amount and percentages
    data_bar=SA_op_stacked_delta; data_scatter=SA_op_stacked_deltanum;
    selidx=1:length(cl_varnames);
    figure;
    subplot(2,1,1)

    b1=bar(categorical(cl_varnames(selidx),cl_varnames(selidx)),data_bar(selidx, :),'stacked');
    grid minor
    ylabel('\Delta TWh/yr')

    yyaxis right
    scatter(categorical(cl_varnames(selidx),cl_varnames(selidx)),data_scatter(:,1),50,cl_RP*.6,'d','filled')
    hold on
    scatter(categorical(cl_varnames(selidx),cl_varnames(selidx)),data_scatter(:,1)+data_scatter(:,2),50,cl_DP*.6,'d','filled') %manually stack scatter
    % legend(t.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
    ylim([-800 750])

    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end
    alpha(b1,0.7)
    ax = gca;
    ax.YAxis(2).Color = 'k';
    ylabel('\Delta # of plants')

    title("Delta change ENERGY for SAother params")

    % SA_Sust other: basic bar plot of delta change PERCENTAGE
    data_bar=SA_op_stacked_deltaprct; data_scatter=SA_op_stacked_deltaprctnum;
    selidx=1:length(cl_varnames);
    subplot(2,1,2)

    b1=bar(categorical(cl_varnames(selidx),cl_varnames(selidx)),data_bar(selidx, :),'stacked');
    ylabel('Percentage change in TWh/yr')
    grid minor
    yyaxis right
    scatter(categorical(cl_varnames(selidx),cl_varnames(selidx)),data_scatter(:,1),50,cl_RP*.6,'d','filled')
    hold on
    scatter(categorical(cl_varnames(selidx),cl_varnames(selidx)),data_scatter(:,1)+data_scatter(:,2),50,cl_DP*.6,'d','filled') %manually stack scatter
    % legend(t.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
    ylim([-146 110])

    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end
    alpha(b1,0.7)
    ax = gca;
    ax.YAxis(2).Color = 'k';
    ylabel('Percentage change in # of plants')
    title("Delta change PERCENTAGE for SAother params")



    %% SA_Sust other: Plot Rotated bar for DELTA Amount SA_Sust other params
    % Cannot make this in one fig, it has to be two figs :E
    % Based on https://nl.mathworks.com/help/matlab/creating_plots/graph-with-multiple-x-axes-and-y-axes.html#responsive_offcanvas
    % Could not replicate this for two plots in one fig case
    data_bar=SA_op_stacked_delta; data_scatter=SA_op_stacked_deltanum;
    runidx= {1:7*3; 7*3+1:7*3*2};
    runtitles=["(+)", "(-)"];
    for i=1:2
        figure
        mytile=tiledlayout(1,1);
        selidx=runidx{i};
        ax1=axes(mytile);
        ax2=axes(mytile);
        b1=barh(ax1,categorical(cl_varnames(selidx),cl_varnames(selidx)),data_bar(selidx, :),'stacked');
        grid minor
        for k = 1:2
            b1(k).FaceColor = mycolors{k};
            b1(k).EdgeColor = mycolors{k};
            b2(k).FaceColor = mycolors{k};
            b2(k).EdgeColor = mycolors{k};
        end
        %
        scatter(data_scatter(selidx,1),categorical(cl_varnames(selidx),cl_varnames(selidx)), 50,cl_RP*.6,'d','filled')
        hold on
        scatter(data_scatter(selidx,1)+data_scatter(selidx,2),categorical(cl_varnames(selidx),cl_varnames(selidx)), 50,cl_DP*.6,'d','filled') %manually stack scatter
        ax2.XAxisLocation = 'top';
        ax2.YAxisLocation = 'right';
        ax2.Color = 'none';
        ax1.Box = 'off';
        ax2.Box = 'off';
        grid minor
        set(ax2,'YTickLabels',[],'XLim',[-80 110])
        title(runtitles(i))
    end
end
%% Write totals to excel and matfile
if writeSummary2file
    save(matfile_summary,'tot_allscen', 'SA_Others_delta','mainidx_end','SA_Sust_idx_end',...
        'SA_Tech_idx_end', 'SA_Fin_idx_end','selmixsust')

    % Write the all scenarios compiled data
    tmp=sprintf("Last updated on: %s",datestr(datetime('now')));
    tot_allscen(size(tot_allscen,1)+1,1)={tmp{:}};
    writetable(tot_allscen,xlsfile_summary,'Sheet',"Compiled");

    % Write the % of plants costing <=0.10
    all_prctfin_num(size(all_prctfin_num,1)+1,1)={tmp{:}};
    writetable(all_prctfin_num,xlsfile_summary,'Sheet',"Prct_fin");


    % Write the SA_Others_delta change data
    SA_Others_delta(size(SA_Others_delta,1)+1,1)={tmp{:}};
    writetable(SA_Others_delta,xlsfile_summary,'Sheet',"SA_Others_delta_R2");
    disp("Written all data to mat and xls")

end