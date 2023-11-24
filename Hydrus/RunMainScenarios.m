%% Run the main scenarios used in the paper
%Run the 3x2x2 energy + 2 geohazrep + 3 geohaz type scenarios using parfor
%Run the 3 Tech2Sust iterations using spmd

clearvars
close all
clc
addpath(genpath(fullfile(pwd,"Hydrus")))
%runname= runHydrus(runname_prefix, rootpath,...
%    scenario, policytype, runsustain, geohazard_scentype, geohazard_select,...
%    makechanges2main, makechanges2cost) % Takes upto 9 PARAMTERS
% runname =<prefix>_<scenariotype>_<policytype>_<constraintype>_<hazardreptype>

% Default nbasin = 101 which has viz = only existing
% manually change nbasin = 102 for vis = exinumScensst+uc
batchRun=1; %If 0 then scenarios are run in for loop else they are run in batches in spmd
trialRun=0;
ncores=feature('numcores');
maxcores=10;
runnameidx=103; %num in filenames
%% Run one scenario - trial works w and wo geohazard
if trialRun
    r_prefix=sprintf('R%d',runnameidx);
    scenario='Full';        % Full or Remain';  %
    policytype=0;           % Policy scenarios 0=trial w mixed, 1= trial w large, 2= large focus, 3= medium focus, 4= mixed focus
    runsustain=0;           % 0 or 1 to turn sustainable constrains off/on
    geohazard_scentype=2;   % For sustainable runs, 0=no geohazard consideration, 1=Costbased, 2= Risk averse, 3= Multi hazard
    geohazard_select=[1 1 1 1]; %Apply EQ-PGA, EQ-Thrust, LS, GLOF hazard risk data
    runHydrus(r_prefix, pwd,...
        scenario, policytype, runsustain, geohazard_scentype, geohazard_select)
    %     runHydrus(r_prefix, pwd,...
    %         scenario, policytype, runsustain, geohazard_scentype, geohazard_select,...
    %         '','interest    = 50/100; ')
    
    %% Try parfor - tested for trial and mixed
    fprintf('\n\n\nTry parfor\n')
    r_prefix='R100';
    scenario='Full';        % Full or Remain';  %
    policytype=4;           % Policy scenarios 0=trial w mixed, 1= trial w large, 2= large focus, 3= medium focus, 4= mixed focus
    runsustain=1;           % 0 or 1 to turn sustainable constrains off/on
    geohazard_scentype=2;   % For sustainable runs, 0=no geohazard consideration, 1=Costbased, 2= Risk averse, 3= Multi hazard
    geohazard_select=[1 1 1 1]; %Apply EQ-PGA, EQ-Thrust, LS, GLOF hazard risk data
    
    numScens=3;runcnt=1;
    pool=parpool;
    parfor (i=1:numScens,min(numScens,ncores))
        runHydrus(strcat(r_prefix,num2str(i)), pwd,...
            scenario, policytype, runsustain, geohazard_scentype, geohazard_select)
    end
    fprintf('\n\n\nEnd Try parfor\n')
end

%% Run one mixed scenario
r_prefix=sprintf('R%d',runnameidx);
   scenario='Full';        % Full or Remain';  %
policytype=2;           % Policy scenarios 0=trial w mixed, 1= trial w large, 2= large focus, 3= medium focus, 4= mixed focus
runsustain=1;           % 0 or 1 to turn sustainable constrains off/on
geohazard_scentype=2; %only run Risk Averse case
geohazard_select=[1 1 1 1]; %apply all 4 disaster data

    runHydrus(r_prefix, pwd,...
        scenario, policytype, runsustain, geohazard_scentype, geohazard_select)

%% New setup Full_Mixed_Tech/Fin to Sust iterations
iterations =["waterconsumption   =1; ";   %Apply water consumption
    "slackflow_constraint=1; ";           %Enable slackflow of x% to ensure natural flow of the river (see slackflow -setting)
    "protarea_constr     =1; " ];         %Apply no dams in protected areas; %
r_prefix=strcat(sprintf('R%d_Tech2Sust',runnameidx), {'+wc';'+wc+eflow';'+wc+eflow+PA'});

scenario={'Full'};        % Full or Remain';  %
policytype=4;           % Policy scenarios 0=trial w mixed, 1= trial w large, 2= large focus, 3= medium focus, 4= mixed focus
runsustain=0;           % 0 or 1 to turn sustainable constrains off/on
geohazard_scentype=2; %only run Risk Averse case
geohazard_select=[1 1 1 1]; %apply all 4 disaster data

runcnt=1;
for i=1:numel(iterations)
    runSetup_wChanges{runcnt}=[r_prefix{i}, pwd, scenario, policytype, runsustain, geohazard_scentype, geohazard_select, {strjoin(iterations(1:i))},' '];
    if ~batchRun
        saveRname_wChanges{runcnt}=runHydrus(runSetup_wChanges{runcnt}{1}, runSetup_wChanges{runcnt}{2},...
            runSetup_wChanges{runcnt}{3}, runSetup_wChanges{runcnt}{4}, runSetup_wChanges{runcnt}{5},...
            runSetup_wChanges{runcnt}{6}, runSetup_wChanges{runcnt}{7}, runSetup_wChanges{runcnt}{8},...
            runSetup_wChanges{runcnt}{9});
    end
    runcnt=runcnt+1;
end

disp("Setup 3 Tech to Sust iterations")


%% New parfor run for Tech to Sust iterative runs
% Loop through constraints
if batchRun
    numScens=numel(runSetup_wChanges);
    CompletedRunxx=cell(numScens,1);
    nWorkers=min([numScens,ncores,maxcores]);

    delete(gcp('nocreate'))
    parpool('local',nWorkers);
    %spmd;   xx=labindex;
     parfor xx=1:numScens;
        CompletedRunxx{xx}=runHydrus(runSetup_wChanges{xx}{1}, runSetup_wChanges{xx}{2},...
            runSetup_wChanges{xx}{3}, runSetup_wChanges{xx}{4}, runSetup_wChanges{xx}{5}, runSetup_wChanges{xx}{6}, runSetup_wChanges{xx}{7}, runSetup_wChanges{xx}{8},runSetup_wChanges{xx}{9});
    end
    delete(gcp)
end
disp("Ran 3 Tech to Sust iterations")

%% New Setup Full/Remain + 3 energy scn + TECH/ECON/SUST_RiskAverse
r_prefix=sprintf('R%d_Energy',runnameidx) ;
geohazard_scentype=2; %only run Risk Averse case
geohazard_select=[1 1 1 1]; %apply all 4 disaster data
runcnt=1;
for scntyp={'Full', 'Remain'}
    for poltyp=2:4
        for susttyp=0:1
            runSetup_Energy{runcnt}=[r_prefix, pwd, scntyp, poltyp, susttyp, geohazard_scentype, geohazard_select];
            if ~batchRun
                saveRname{runcnt}=runHydrus(runSetup_Energy{runcnt}{1}, runSetup_Energy{runcnt}{2},...
                    runSetup_Energy{runcnt}{3}, runSetup_Energy{runcnt}{4}, runSetup_Energy{runcnt}{5},...
                    runSetup_Energy{runcnt}{6}, runSetup_Energy{runcnt}{7});
            end
            runcnt=runcnt+1;
        end
    end
end
disp("Setup 3x2x2 energy scenarios")


%% New Setup Full_Mixed_Sust + 2 Geohazard Rep iterations
r_prefix=sprintf('R%d_HazRep',runnameidx) ;
runcnt=1;
scenario={'Full'};        % Full or Remain';  %
policytype=4;           % Policy scenarios 0=trial w mixed, 1= trial w large, 2= large focus, 3= medium focus, 4= mixed focus
runsustain=1;           % 0 or 1 to turn sustainable constrains off/on
geohazard_select=[1 1 1 1]; %apply all 4 disaster data

for geohazscntyp=[1,3]
    runSetup_GeoHaz{runcnt}=[r_prefix, pwd, scenario, policytype, runsustain, geohazscntyp, geohazard_select];
    if ~batchRun
        saveRname{runcnt}=runHydrus(runSetup_GeoHaz{runcnt}{1}, runSetup_GeoHaz{runcnt}{2},...
            runSetup_GeoHaz{runcnt}{3}, runSetup_GeoHaz{runcnt}{4}, runSetup_GeoHaz{runcnt}{5},...
            runSetup_GeoHaz{runcnt}{6}, runSetup_GeoHaz{runcnt}{7});
    end
    runcnt=runcnt+1;
end
disp("Setup 2 geohazard rep variations")

%% New Setup Full_Mixed_Sust_RiskAverse + One Geohazard at a time iterations
r_prefix=sprintf('R%d_HazTyp',runnameidx) ;
geohazard_types=  [1 1 0 0
    0 0 1 0
    0 0 0 1];   %apply all 4 disaster data
scenario={'Full'};          % Full or Remain';  %
policytype=4;               % Policy scenarios 0=trial w mixed, 1= trial w large, 2= large focus, 3= medium focus, 4= mixed focus
runsustain=1;               % 0 or 1 to turn sustainable constrains off/on
geohazard_scentype=2;       %only run Risk Averse case

for i=1:size(geohazard_types,1)
    runSetup_GeoHaz{runcnt}=[r_prefix, pwd, scenario, policytype, runsustain, geohazard_scentype, geohazard_types(i,:)];
    if ~batchRun
        saveRname{runcnt}=runHydrus(runSetup_GeoHaz{runcnt}{1}, runSetup_GeoHaz{runcnt}{2},...
            runSetup_GeoHaz{runcnt}{3}, runSetup_GeoHaz{runcnt}{4}, runSetup_GeoHaz{runcnt}{5},...
            runSetup_GeoHaz{runcnt}{6}, runSetup_GeoHaz{runcnt}{7});
    end
    runcnt=runcnt+1;
end
disp("Setup 3 variations of one geohazard at a time")
    
%% New Parfor run for main scens: 3 energy (12 total), + 2 geohazard reps, + 3 geohaz type
% Loop through constraints
%runSetup= runSetup_Energy; %ncores=8; % took 4 hrs on the new hpc
runSetup= runSetup_GeoHaz; %ncores=5;

numScens=numel(runSetup);
nWorkers=min([numScens,ncores,maxcores]);
CompletedRun=strings(numScens,1);
%selRun=[3:6]; 15
selRun= 1:numScens; %[9:17];
fprintf("Running %d scenarios in parfor mode\n", numel(selRun))
if batchRun
    delete(gcp('nocreate'))
    parpool('local', nWorkers);
    %spmd;   xx=labindex;
    parfor (xx=selRun,nWorkers)
               %for xx=1:numScens;
         CompletedRun(xx)=runHydrus(runSetup{xx}{1}, runSetup{xx}{2},...
            runSetup{xx}{3}, runSetup{xx}{4}, runSetup{xx}{5}, runSetup{xx}{6}, runSetup{xx}{7});
    end
    delete(gcp)
end
disp("Ran the main scenarios")
