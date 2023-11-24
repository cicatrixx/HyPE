%% Aggregate monthly TS from SimulateHPproduction,to annual TS and performance stats and making exploratory plots
% There are 50 portfolios (1 hist+3ssp*2tf*4corner=25 * 2 for Tech/sust)
% run under 25Qs
% Get basin totals
% Files are saved separately for tech and sust pfs in GWh

% For each pf, qf combo calculate the annual perf statistics for the
% simulation years (hist, mid fut, far fut)
clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

%fnames w/o extension as it is added later
ifmatHPgenTS=fullfile(rootof,'FutScenarios_HPgen_allTS03');

ofmatHPgenTS_basintot=fullfile(rootof,'FutScenarios_HPgen_basintot03'); 
ofmatHPgenTS_prj=fullfile(rootof,'FutScenarios_HPgen_prjstat03');

%% Get yearlyTS and monthlyavgTS for each PF and Q combo
logf = fullfile(rootof, sprintf('Run compile HP performance on %s.log',string(datetime('now','Format','dd-MMM-yyyy'))));
diary(logf)
tic

for pot=[1,4]
    % Load HP performance TS
    currpotname=  newpots4{pot}; % pottypes2{pot};
    load(strcat(ifmatHPgenTS,"_",currpotname,".mat"))
    % Chnage var loading for -Fin type pot
if rem(pot,2)==0
    PactualTS_GWh_4histQ=PactualTS_GWh_4histQ_fin;
    PactualTS_GWh_4futQ=PactualTS_GWh_4futQ_fin;
    clear PactualTS_GWh_4histQ_fin PactualTS_GWh_4futQ_fin
end

    nqfs=length(cl_Qnames);
    cl_pfnames=strcat(extractBefore(pfnames(:,1),"_Tech"),"_",currpotname);
    for qf=1:nqfs %13 = all Qs in near fut
        npfs=length(PactualTS_GWh_4histQ);
        fprintf("Compiling HPgeneration TS for %s pot %s #%d out of %d\n",currpotname,cl_Qnames{qf},qf,nqfs)

         % sel right tf for fut - mid or far
        if qf<=13 %MF
            sel_fut_ts=futmon_ts1;
        else %FF
            sel_fut_ts=futmon_ts2;
        end

        for pf=1:npfs
            if qf==1 % for historical Q
                % Get mon and annual TS at project and basin level
                mytimestamp=histmon_ts;
                tmp_monTS=splitvars(timetable(mytimestamp, PactualTS_GWh_4histQ{pf}(:,:,qf)','VariableNames',"Projects"));
                tmp_monTS.basintot=sum(tmp_monTS{:,:},2);
                tmp_annualTS=table2array(retime(tmp_monTS, 'yearly', 'sum'));  % This contains both project and basintot
                hist_annualTS{pf}(:,:,qf)=tmp_annualTS; % save this as useful for sub-basin analysis

                % Save basin level separately for easier retrival
                hist_monTS_basin(pf,:,qf)=tmp_monTS.basintot;
                hist_annualTS_basin(pf,:,qf)=tmp_annualTS(:,end);
            else
                % Get mon and annual TS at project and basin level
                mytimestamp=sel_fut_ts;
                tmp_monTS=splitvars(timetable(mytimestamp, PactualTS_GWh_4futQ{pf}(:,:,qf)','VariableNames',"Projects"));
                tmp_monTS.basintot=sum(tmp_monTS{:,:},2);
                tmp_annualTS=table2array(retime(tmp_monTS, 'yearly', 'sum'));  % This contains both project and basintot
                fut_annualTS{pf}(:,:,qf)=tmp_annualTS; % save this as useful for sub-basin analysis

                % Save basin level separately for easier retrival
                fut_monTS_basin(pf,:,qf)=tmp_monTS.basintot;
                fut_annualTS_basin(pf,:,qf)=tmp_annualTS(:,end);
            end
            % Basin level outputs that are same length for hist and fut so can be put in one table
            % Get LT monthly stats - only do this at basin scale for now
            % it does not make sense to calc any medians in a project level, might need sub-basin level but that will have to be done separately anyways
            tmp_LTmon=groupsummary(tmp_monTS,"mytimestamp","monthname",["mean","median","min","max","std","var"],"basintot");

            % These do not need to be separated for hist and fut because
            % they have the same time frame.
            % Archive table outputs into a large structure
            LTmon_basin.mean(pf,:,qf) = tmp_LTmon.mean_basintot';
            LTmon_basin.median(pf,:,qf) = tmp_LTmon.median_basintot';
            LTmon_basin.min(pf,:,qf) = tmp_LTmon.min_basintot';
            LTmon_basin.max(pf,:,qf) = tmp_LTmon.max_basintot';
            LTmon_basin.std(pf,:,qf) = tmp_LTmon.std_basintot';
            LTmon_basin.var(pf,:,qf) = tmp_LTmon.var_basintot';
            
            % Get LT annual stats, basin TS is at the last column of LTannual_basin
            LTannual_basin.mean(pf,1,qf)=mean(tmp_annualTS(:,end));
            LTannual_basin.median(pf,1,qf)=median(tmp_annualTS(:,end));
            LTannual_basin.max(pf,1,qf)=max(tmp_annualTS(:,end));
            LTannual_basin.min(pf,1,qf)=min(tmp_annualTS(:,end));
            LTannual_basin.std(pf,1,qf)=std(tmp_annualTS(:,end));
            LTannual_basin.var(pf,1,qf)=var(tmp_annualTS(:,end));


             %% Add annual percentile calculations
            if qf==1 % for historical Q
                basin_annualTS=hist_annualTS_basin(pf,:,qf);
            else
                basin_annualTS= fut_annualTS_basin(pf,:,qf);
            end
            LTannual_basin.prctile10(pf,qf)=findPrctile(basin_annualTS,10);
            LTannual_basin.prctile90(pf,qf)=findPrctile(basin_annualTS,90);
          %  fprintf("EvaluatedTS stats for %s potential Q#%d PF#%d\n",currpotname,qf,pf)
        end
    end
    %% Save files
    % in these files, hist_/fut_ indicates data for histQ and futQs. data is
    % then stored in matrices for each portfolio for each futQs
    save(strcat(ofmatHPgenTS_basintot,"_",currpotname,".mat"),'LTannual_basin',"LTmon_basin",...
        "hist_annualTS_basin", "hist_monTS_basin","fut_annualTS_basin", "fut_monTS_basin","cl_pfnames","cl_Qnames")
    save(strcat(ofmatHPgenTS_prj,"_",currpotname,".mat"),...
        "hist_annualTS", "fut_annualTS","cl_pfnames","cl_Qnames")

    fprintf("Saved basin level HP TS stats for %s potential\n", currpotname)
    clear tmp_* LT* *_annualTS* *_monTS* PactualTS*
end
disp("Completed evaluation of HP TS stats")
diary off
toc
%             %% Add annual percentile calculations 
% 
% for pot=1
%     load(strcat(ofmatHPgenTS_basintot,"_",pottypes2{pot},".mat"),'LTannual_basin',...
%         "hist_annualTS_basin","fut_annualTS_basin")
%     npfs=size(fut_annualTS_basin,1);
%     nqfs=size(fut_annualTS_basin,3);
%     for qf=1:nqfs %13 = all Qs in near fut
%         %fprintf("Compiling HPgeneration TS for Q #%d out of %d\n",qf,nqfs)
%         for pf=1:npfs
%             if qf==1 % for historical Q
%                 basin_annualTS=hist_annualTS_basin(pf,:,qf);
%             else
%                 basin_annualTS= fut_annualTS_basin(pf,:,qf);
%             end
%             LTannual_basin.prctile10(pf,qf)=findPrctile(basin_annualTS,10);
%             LTannual_basin.prctile90(pf,qf)=findPrctile(basin_annualTS,90);
%             %matprctile10(pf,qf)=prctile(basin_annualTS,90); % decided to avoid the interpolation done here            
%         end
%     end
% end
% disp("Evaluated percentiles")
% save(strcat(ofmatHPgenTS_basintot,"_",pottypes2{pot},".mat"),'LTannual_basin','-append')
