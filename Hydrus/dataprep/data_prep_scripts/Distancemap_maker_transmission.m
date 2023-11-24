%% Distance map 
% Distance to nearest transmission line (# of cells)
% Voor NAM duurt deze run 3 dagen
% Original by David Gernaat
% SD changed:
% - improved buffer window creation and indexing
% - converted to function

% Run time 2hrs w parfor
clear all
%close all

cellsz=500; % in m
search_rad=500*1e3; %in m %max in data was 479 so cant reduce this for future run
win=search_rad/cellsz; % # of cells

%% Load DIR and Powerlines.mat for selected continent
disp('Loading dataset to assess distance to');
load(fullfile(pwd,'data','UI','data','Powerlines.mat'), 'data');

disp('Loading outside bounds');
load(fullfile(pwd,'data','UI','data','Basin','Basin_551.mat'),'outside')
data(outside)=0;
Infeature=data;


Infeature=logical(Infeature);
%% Create subset for test run w ASIA file
rl=1;
testc=1300:1800;
testr=400:900;
if rl==2
    Infeature=Infeature(testr, testc);
    outside=logical(outside(testr, testc));
end
figure(1);clf;imagesc(Infeature);axis image
%figure(2);clf;imagesc(outside);axis image
%%
Dis=getDistanceMap(Infeature,outside,win);

%overlay feature into map
[r, c]=find(Infeature);
hold all
scatter(c,r,'.k')
legend('Transmission line')
c=colorbar;
c.Label.String='Distance to nearest transmission line (# of cells)';
%% Save
disp('Save')
data=Dis;
matfile = fullfile(pwd,'data','UI','data','Dis2Transmission.mat');
save(matfile,'data');