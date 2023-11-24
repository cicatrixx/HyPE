function adir = reversePatrickDirs(fdir)

% Load directional metrics for Patrick dirs
run('defdirs.m' )

fprintf('Reversing flow directions\n');
adir = zeros(size(fdir), 'int8');
for d=1:9
    idx = fdir==d;
    adir(idx) = antidir(d);
end;

%% Test
% 
% fname='G:/PaperData/DavidGernaat/Sanita_model_package/data/ASIA/Basin/Basin_5.mat';
% load(fname,'fdir','adir')
% xx=reversePatrickDirs(fdir);
% all(xx(:)==adir(:))
