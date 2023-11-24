%% Run one scenario - trial works w or wo geohazard
% Call script from HPmodel folder
clearvars
close all
clc
addpath(genpath(fullfile(pwd,"Hydrus")))

runnameidx=103;
r_prefix=sprintf('TRIAL_Fut_HPC %d',runnameidx);
scenario='Full';        % Full or Remain';  %
policytype=0;           % Policy scenarios 0=trial w mixed, 1= trial w large, 2= large focus, 3= medium focus, 4= mixed focus
runsustain=0;           % 0 or 1 to turn sustainable constrains off/on
geohazard_scentype=2;   % For sustainable runs, 0=no geohazard consideration, 1=Costbased, 2= Risk averse, 3= Multi hazard
geohazard_select=[1 1 1 1]; %Apply EQ-PGA, EQ-Thrust, LS, GLOF hazard risk data
runHydrus(r_prefix, pwd,...
    scenario, policytype, runsustain, geohazard_scentype, geohazard_select,'','',...
    fullfile(pwd,'data','ASIA','Basin_UIB','Fut_designQs','ssp245_AWI-CM-1-1-MR_LTavgs_2036_2065.mat'))
