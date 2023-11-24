
function Dis=getDistanceMap(Idata,Outside,win)
% Adapted from David Gernaat et al 2017
% assigns minimum distance of half a cell. If no feature found within
% search radius, then dis=dsmax*1e3

%% Add buffers for searchwindow
[nr,nc]=size(Idata);
%Create extra window for overlaying w repmat at the borders
Dis_win = zeros((nr+2*win),(nc+2*win),'single');
Idata_win = zeros((nr+2*win),(nc+2*win),'int8');
Idata_win(win+1:nr+win,win+1:nc+win)=logical(Idata);

Outside_win = ones((nr+2*win),(nc+2*win),'int16');
Outside_win(win+1:nr+win,win+1:nc+win)=logical(Outside);
disp('Create rasters with borders')

%% Calculate distance
disp('Calculate distance')
pold=0;
N=nr*nc;
i=0;

% Build repmat cell distance map
dx = repmat(abs(-win:win), [2*win+1, 1]);
dy = dx';
distance_win = sqrt(dx.^2+dy.^2);
dsmax = max(distance_win,[],'all');
cvals=(win+1):(nc+win);

tic
parfor r=(win+1):(nr+win)
    for c=cvals; %(win+1):(nc+win)
        %NaNs for sea
        if Outside_win(r,c)==1; Dis_win(r,c)=0; continue; end
        
        % Build window around target cell
        firstrow = r-win;
        lastrow  = r+win;
        firstcol = c-win;
        lastcol  = c+win;
        
        Pwin = Idata_win(firstrow:lastrow,firstcol:lastcol);
        
        %figure(3);clf;imagesc(log(demwin));axis image;
        
        %Calc distance using repmat dis table
        near_w = distance_win(Pwin==1);
        dis = max(0.5,min(near_w,[],'all'));
        if isempty(dis)==true; dis=dsmax*1e3; end
        Dis_win(r,c)=dis;
        % Report
        %         i=i+1;
        %         p = fix(1000*i/N);
        %         if p>pold
        %             fprintf('#%d Progress %d%%.\n',rl,p);
        %             pold=p;
        %         end;
    end
end
toc

disp('Distance map created!')

%% Cut correct window for output
Dis = Dis_win(win+1:(nr+win),win+1:(nc+win)); % in cells
figure;clf;imagesc(Dis);axis image;colormap(jet);

end