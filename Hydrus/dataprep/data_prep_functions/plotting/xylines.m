function xylines(pts,LineSpec,X_Y)
% Helper function to plot multiple xlines or ylines
X_Y=upper(X_Y);
if X_Y=="X"
    for i=1:numel(pts)
        xline(pts(i),LineSpec);
    end
elseif X_Y=="Y"
    for i=1:numel(pts)
        yline(pts(i),LineSpec);
        
    end
else
    disp("!!! Not clear if x or y line !!!")
end