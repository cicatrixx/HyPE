%function demwin = getHighResDEMwin(, winsize100m)
winsize100m=30;
winsize500m=winsize100m*100/500;

%% Get 500 win
%create window in low-res dem
[nr500,nc500]=size(Z);
firstrow = max(1,ro(k)-winsize500m);
lastrow  = min(nr500,ro(k)+winsize500m);
firstcol = max(1,co(k)-winsize500m);
lastcol  = min(nc500,co(k)+winsize500m);
Zwinlow=Z(firstrow:lastrow,firstcol:lastcol);

%% Get high-resolution window data
load('G:\SurfDrive\GitConnect\data\UI\data\Z100m.mat')
Z_100m=data;
Rw_100m=Rw100;
[nr100,nc100]=size(Z_100m);
%get RC from XY coordinate for 100m DEM
[r_dis,c_dis] = worldToDiscrete(Rw_100m,xOut(k),yOut(k));

%create window in high-res dem
firstrow = max(1,r_dis-winsize100m);
lastrow  = min(nr100,r_dis+winsize100m);
firstcol = max(1,c_dis-winsize100m);
lastcol  = min(nc100,c_dis+winsize100m);
Zwinhigh = Z_100m(firstrow:lastrow,firstcol:lastcol);

%% Plot both resolution data and r,c point to check 
figure
subplot(2,1,1)
imagesc(Z)
hold on 
plot(co(k),ro(k),'.r')
axis image

subplot(2,1,2)
imagesc(Z_100m)
hold on 
plot(c_dis,r_dis,'.r')
axis image

%% Plot Z windows to check
figure
subplot(2,1,1)
imagesc(Zwinhigh)
axis image

subplot(2,1,2)
imagesc(Zwinlow)
axis image
