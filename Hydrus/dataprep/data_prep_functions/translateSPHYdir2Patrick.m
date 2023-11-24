function fdir = translateSPHYdir2Patrick(dir)
% Convert ESRI flow direction definition to Patrick direction metrics.
% Similar to set up in David Gernaat's codes in BIL_conv.m

fprintf('Translating flow directions\n');
%From def dirs ESRI to Patrick dirs mapping
mapdir = [
    %Patrick    %ESRI     %SPHY %Nominal
    1           64        8     %NORTH
    2           16        4     %WEST
    3           1         6     %EAST
    4           4         2     %SOUTH
    5           32        7     %NWEST
    6           128       9     %NEAST
    7           8         1     %SWEST
    8           2         3     %SEAST
    9           247       5     %SELF/SINK % also 255 for ESRI
    ];

SPHYdirs = mapdir(:,3);
Patdirs = mapdir(:,1);

fdir = zeros(size(dir), 'int8');
% Fill Patrick codes 1-9 in cells w equivalent SPHYdirs
for d=1:9
    idx = dir(:) == SPHYdirs(d);
    fdir(idx) = d;
end;

%% CHeck: Compare SPHY and Patrick dirs
checkDirs=1;
if checkDirs
    % SPHYdirs unique counts
    value_counts_dir = countUniques(dir);
    
    % Patrickdirs unique counts
    value_counts_fdir = countUniques(fdir);
    
    cmpcnts=table(value_counts_dir(SPHYdirs,:), value_counts_fdir,'VariableNames', {'SPHY', 'Patrick'});
    if all(value_counts_dir(SPHYdirs,2) == value_counts_fdir(:,2))
        disp("Count in direction categories match")
    end
end