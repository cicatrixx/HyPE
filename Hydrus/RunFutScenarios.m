%% Run all future scenarios
clear
close all
clc
addpath(genpath(fullfile(pwd,"Hydrus")))
batchRun=1; %If 0 then scenarios are run in for loop else they are run in parallel batches
runnameidx=103; %num in filenames
r_prefix=sprintf('FutR%d',runnameidx);
scenario='Full';            % Full or Remain';  %
policytype=4;               % Policy scenarios 0=trial w mixed, 1= trial w large, 2= large focus, 3= medium focus, 4= mixed focus
geohazard_scentype=2;       % only run Risk Averse case
geohazard_select=[1 1 1 1]; % apply all 4 disaster data
loopnm={'Tech-Fin', 'Sust'};
runtype=1;                  % 2=sust only or 1=tech only

%Qchange="sel_design_flows runHydrus= {'%s', '%s', '%s'};";  %RP_Qs{ii},largeDP_Qs{jj},smallDP_Qs
ncores=feature('numcores');
%% Get paths to future files
Qpath= fullfile(pwd,'data','ASIA','Basin_UIB','Fut_designQs');
[fpaths_Q500m, fname]=getPaths2Files(Qpath,"ssp*");  
fname=strrep(fname,'.mat','');
numScens=length(fpaths_Q500m);

%% New PARFOR run for future SUSTAIN scenarios: 24 cases (3RCPs, 4Models, 2TFs) x 2 (Tech-Fin and Sust-Fin)
% Get num of parallel workers as the minimum of num of cores, num of scenarios or max cores to use to ensure parallelization does not exceed memory
maxcores=3; % for ram-128GB HPC, can run about 8 scenarios in parallel as i need max 16GB for one scen % the tech_mix scen are is the heaviest
% did not work w 7 
numScens=1;
nWorkers=min([numScens,ncores,maxcores]);
CompletedRun_sust=strings(numScens,1);
selRun=[14]; %:numScens; %8 remains
if runtype==2
    if batchRun
        delete(gcp('nocreate'))
        parpool('local', nWorkers);
        fprintf("Running %d scenarios in parfor mode\n", numel(selRun))
        for runsustain=1           % 0 or 1 to turn sustainable constrains off/on
            %cidx=runsustain*numScens;
            parfor (xxy=1,nWorkers)
                xx=selRun(xxy);
                %runname= runHydrus(runname_prefix, rootpath,...
                %    scenario, policytype, runsustain, geohazard_scentype, geohazard_select,...
                %    makechanges2main, makechanges2cost, path2Q) % Takes upto 10 PARAMTERS
                % runname =<prefix>_<scenariotype>_<policytype>_<constraintype>_<hazardreptype>

                runname_prefix=[r_prefix,'_', fname{xx}];
                CompletedRun_sust(xxy)=runHydrus(runname_prefix,...
                    pwd, scenario, policytype, runsustain, geohazard_scentype, geohazard_select,'','',...
                    fpaths_Q500m{xx});
            end
            disp(CompletedRun_sust)
            delete(gcp)
        end

        fprintf("Ran the Future %s scenarios", loopnm{runsustain+1})
    end
end

disp(">>>>>>>>>>>>>>>>>>>>>>>>>EOF>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

%% New PARFOR run for future TECH scenarios: 24 cases (3RCPs, 4Models, 2TFs) x 2 (Tech-Fin and Sust-Fin)
% Get num of parallel workers as the minimum of num of cores, num of scenarios or max cores to use to ensure parallelization does not exceed memory
maxcores=2; % for ram-128GB HPC, can run about 8 scenarios in parallel as i need max 16GB for one scen % the tech_mix scen are is the heaviest
% did not work w 7 or 4 , perhaps the tech-fin scens are too heavy

nWorkers=1; %min([numScens,ncores,maxcores]);
CompletedRun_tech=strings(numScens,1);
selRun=24;%numScens; %8 remains
if runtype==1
    if batchRun
        delete(gcp('nocreate'))
        parpool('local', nWorkers);
        fprintf("Running %d scenarios in parfor mode\n", numel(selRun))
        for runsustain=0           % 0 or 1 to turn sustainable constrains off/on
            parfor (xxy=1, nWorkers)
                xx=selRun(xxy);

                %runname= runHydrus(runname_prefix, rootpath,...
                %    scenario, policytype, runsustain, geohazard_scentype, geohazard_select,...
                %    makechanges2main, makechanges2cost, path2Q) % Takes upto 10 PARAMTERS
                % runname =<prefix>_<scenariotype>_<policytype>_<constraintype>_<hazardreptype>

                runname_prefix=[r_prefix,'_', fname{xx}];
                CompletedRun_tech(xxy)=runHydrus(runname_prefix,...
                    pwd, scenario, policytype, runsustain, geohazard_scentype, geohazard_select,'','',...
                    fpaths_Q500m{xx});
            end
            disp(CompletedRun_tech)
            delete(gcp)
        end

        fprintf("Ran the Future %s scenarios", loopnm{runsustain+1})
    end
end
disp(">>>>>>>>>>>>>>>>>>>>>>>>>EOF>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

