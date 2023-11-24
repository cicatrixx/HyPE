%% Run the SA scenarios used in the paper
% for ram-128GB HPC, can run about 8 scenarios in parallel as i need max 16GB
% SA_Q for RP3 is intensive needs to be run w fewer cores
% for the tech_mix scen that is the heaviest
clear
close all
clc
addpath(genpath(fullfile(pwd,"Hydrus")))
%runname= runHydrus(runname_prefix, rootpath,...
%    scenario, policytype, runsustain, geohazard_scentype, geohazard_select,...
%    makechanges2main, makechanges2cost) % Takes upto 9 PARAMTERS
% runname =<prefix>_<scenariotype>_<policytype>_<constraintype>_<hazardreptype>

% Default nbasin = 101 which has viz = only existing
% manually change nbasin = 102 for vis = exist+uc
batchRun=1; %If 0 then scenarios are run in for loop else they are run in batches in spmd

runnameidx=103; %num in filenames
r_prefix=sprintf('R%d_SA',runnameidx);
scenario='Full';          % Full or Remain';  %
policytype=4;               % Policy scenarios 0=trial w mixed, 1= trial w large, 2= large focus, 3= medium focus, 4= mixed focus
runsustain=1;               % 0 or 1 to turn sustainable constrains off/on
geohazard_scentype=2;       % only run Risk Averse case
geohazard_select=[1 1 1 1]; % apply all 4 disaster data
runSA_Qs=0;
runSA_Sust=0;
runSA_TechFin=1;

ncores=feature('numcores');
maxcores=7;

%% SA runs for Q for the Full_Mixed_Sust_Risk Averse case
if runSA_Qs
    RP_Qs =      {'Q30', 'Q40', 'Q50'};
    largeDP_Qs = {'Q25', 'Q30', 'Q40'};
    smallDP_Qs = {'Q70', 'Q80', 'Q90'};
    runcnt=1;

    for ii=1:numel(RP_Qs)
        for jj=1:numel(largeDP_Qs)
            for kk=1:numel(smallDP_Qs)
                genprefix={strcat(r_prefix, sprintf('_R%dlDP%dsDP%d',ii,jj,kk))};
                genQ={sprintf("sel_design_flows = {'%s', '%s', '%s'};",RP_Qs{ii},largeDP_Qs{jj},smallDP_Qs{kk})};
                runSetup_Qs{runcnt}=[genprefix, genQ];
                if ~batchRun
                    saveRname_Qs{runcnt}=runHydrustmp(runSetup_Qs{runcnt}{1}, ...
                        pwd, scenario, policytype, runsustain, geohazard_scentype, geohazard_select,...
                        runSetup_Qs{runcnt}{2});
                end
                runcnt=runcnt+1;
            end
        end
    end
    disp("Setup SA_Q runs")

    %% New SPMD run for SA for Q scens: 27 cases
    % % Loop through constraints
    % numScens=numel(runSetup_Qs);
    % ncores=feature('numcores');
    % nWorkers=min(numScens,ncores);
    % CompletedRun=strings(numScens,1);
    % selRun=1:27;
    %
    % if numScens>nWorkers
    %     breakRuns=0:nWorkers:numScens;
    % else
    %     breakRuns=0;
    % end
    %
    % fprintf("Running %d scenarios in parfor mode\n", numel(selRun))
    % if batchRun
    %     delete(gcp('nocreate'))
    %     parpool('local', nWorkers);
    %     for i=1:length(breakRuns)
    %         spmd;   xx=labindex+breakRuns(3);
    %             %parfor (xx=selRun,nWorkers)
    %             CompletedRun(xx)=runHydrustmp(runSetup_Qs{xx}{1},...
    %                 pwd, scenario, policytype, runsustain, geohazard_scentype, geohazard_select,...
    %                 runSetup_Qs{xx}{2});
    %         end
    %     end
    %     delete(gcp)
    % end
    % disp("Ran the SA Q scenarios")

    %% New PARFOR run for SA for Q scens: 27 cases
    % Loop through constraints
    numScens=numel(runSetup_Qs);
    %ncores=6;%feature('numcores');
maxcores=4;
    nWorkers=min([numScens,ncores,maxcores]);
    %CompletedRun=strings(numScens,1); [25:27];%
    selRun=22:27;%

    if batchRun
        delete(gcp('nocreate'))
        parpool('local', nWorkers);
        %for i=1:(length(breakRuns)-1)
        %selRun=(breakRuns(i)+1):breakRuns(i+1);
        fprintf("Running %d scenarios in parfor mode\n", numel(selRun))
        parfor (xx=selRun,nWorkers)
            CompletedRun(xx)=runHydrus(runSetup_Qs{xx}{1},...
                pwd, scenario, policytype, runsustain, geohazard_scentype, geohazard_select,...
                runSetup_Qs{xx}{2},' ');
        end
        CompletedRun
        %end
        delete(gcp)
    end
    disp("Ran the SA Q scenarios")

end

%% Other SA runs: 7x2
%Load excel file w param list
SAruns = readtable(fullfile(pwd,"SensitivityAnalysis_Setup.xlsx"),'Range','2:9'); %
runcnt=1;
SAruns.Min(2:end) =SAruns.Min(2:end)/100;
SAruns.Max(2:end) =SAruns.Max(2:end)/100;

%Create SA eqns
for SAtype={'Min','Max'}
    for parm=1:height(SAruns)
        genprefix={strjoin([r_prefix, SAtype, SAruns.Parameter{parm}],"_")};
        parmchange=sprintf("%s = %.2f;",SAruns.Parameter{parm},SAruns.(SAtype{1})(parm));

        if strcmp(SAruns.ParameterLoc{parm},'MainConfig')
            changemain=parmchange;
            changecost=' ';
        elseif strcmp(SAruns.ParameterLoc{parm},'CostConfig')
            changemain=' ';
            changecost=parmchange;
        end
        runSetup_SAs{runcnt}=[genprefix,  changemain, changecost];
        if ~batchRun
            saveRname_SAs{runcnt}=runHydrustmp(runSetup_SAs{runcnt}{1},...
                pwd, scenario, policytype, runsustain, geohazard_scentype, geohazard_select,...
                runSetup_SAs{runcnt}{2},runSetup_SAs{runcnt}{3});
        end
        runcnt=runcnt+1;
    end
end

%% New PARFOR run for SA for other vars in Sust case: 14 cases
if runSA_Sust
    runsustain=1;
    % Loop through constraints
    runSetup=runSetup_SAs;
    numScens=numel(runSetup);
    %ncores=3;%feature('numcores');
    nWorkers=min(numScens,ncores);
    CompletedRun=strings(numScens,1);
    selRun=1:numScens;

    if batchRun
        delete(gcp('nocreate'))
        parpool('local', nWorkers);
        fprintf("Running %d scenarios in parfor mode\n", numel(selRun))
        parfor (xx=selRun,nWorkers)
            CompletedRun(xx)=runHydrus(runSetup{xx}{1},...
                pwd, scenario, policytype, runsustain, geohazard_scentype, geohazard_select,...
                runSetup{xx}{2},runSetup{xx}{3});
        end
        CompletedRun
        %end
        delete(gcp)
    end

    disp("Ran the SA Sust scenarios")
end
%% New PARFOR run for SA for other vars in Tech-Fin scen: 14 cases
if runSA_TechFin
    runsustain=0;
    % Loop through constraints
    runSetup=runSetup_SAs;
    %{2}=runSetup_SAs{5};
    numScens=numel(runSetup);
    maxcores=7;
    nWorkers=min([numScens,ncores,maxcores]);
    %CompletedRun=strings(numScens,1);
    selRun=1:numScens; %8 remains

    if batchRun
        delete(gcp('nocreate'))
        parpool('local', nWorkers);
        fprintf("Running %d scenarios in parfor mode\n", numel(selRun))
        parfor (xx=selRun,nWorkers)
            CompletedRun(xx)=runHydrus(runSetup{xx}{1},...
                pwd, scenario, policytype, runsustain, geohazard_scentype, geohazard_select,...
                runSetup{xx}{2},runSetup{xx}{3});
        end
        CompletedRun
        %end
        delete(gcp)
    end

    disp("Ran the SA Tech-Fin scenarios")
end