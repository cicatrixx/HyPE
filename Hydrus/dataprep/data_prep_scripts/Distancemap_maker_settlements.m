%% Distance map
% Voor NAM duurt deze run 3 dagen
% Original by David Gernaat
% SD changed:
% - improved buffer window creation and indexing
% Elapsed time is 2094.657246 seconds.

clear all
close all

cellsz=500; % in m
search_rad=300*1e3; %in m %max in data was 120km so can reduce this for future run
win=search_rad/cellsz % # of cells

%% Load DIR and Powerlines.mat for selected continent
rl=	1;

disp('Loading dataset to assess distance to');
load(fullfile(pwd,'data','UI','data','Settlements.mat'), 'data');

disp('Loading outside bounds');
load(fullfile(pwd,'data','UI','data','Basin','Basin_551.mat'),'outside')
data(outside)=0;
Infeature=data;

Infeature=logical(Infeature);
%% Create subset for test run w ASIA file
testc=1300:1800;
testr=400:900;
if rl==2
    Infeature=Infeature(testr, testc);
    outside=logical(outside(testr, testc));
end
figure(1);clf;imagesc(Infeature);axis image
figure(2);clf;imagesc(outside);axis image

%%
Dis=getDistanceMap(Infeature,outside,win);
%overlay feature into map
[r, c]=find(Infeature);
hold all
scatter(c,r,'.k')
legend('Settlements')
c=colorbar;
c.Label.String='Distance to nearest settlement (# of cells)';

%% Save
disp('Save')
data=Dis;
matfile = 'G:\SurfDrive\GitConnect\data\UI\data\Dis2Settlement.mat';
save(matfile,'data');