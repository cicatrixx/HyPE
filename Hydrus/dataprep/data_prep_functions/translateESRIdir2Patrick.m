function fdir = translateESRIdir2Patrick(dir)
% Convert ESRI flow direction definition to Patrick direction metrics.
% Code copied as if from David Gernaat's codes in BIL_conv.m

fprintf('Translating flow directions\n');

% Load directional metrics for Patrick dirs
run('defdirs.m' )
fdir = zeros(size(dir), 'int8');

% Fill Patrick codes 1-8 in cells w equivalent ESRIdirs
for d=1:8
    idx = dir(:) == ESRIdirs(d);
    fdir(idx) = d;
end;

% Fill Patrick code 9 in cells w sea/sink/self points in ESRIdirs
idx_sea = dir(:) == 247; % Sea points
fdir(idx_sea) = 9; % SELF

idx_sink = dir(:) == 255; % Sink points
fdir(idx_sink) = 9; % SELF

idx_outlet = dir(:) == 0; % nan points
fdir(idx_outlet) = 9; % no data


%% CHeck: Compare SPHY and Patrick dirs
checkDirs=1;
if checkDirs
    % ESRI dirs unique counts
    value_counts_dir = countUniques(dir);
    
    % Patrickdirs unique counts
    value_counts_fdir = countUniques(fdir);
    
    cmpcnts=table(value_counts_dir, value_counts_fdir,'VariableNames', {'ESRI', 'Patrick'})
%     if all(value_counts_dir(:,2) == value_counts_fdir(:,2))
%         disp("Count in direction categories match")
%     end
end

end

