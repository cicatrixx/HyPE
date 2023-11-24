%% Pop resizer

clear all

root{1} = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\NAM';
root{2} = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global\EUR';
root{3} = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global\SAM';
root{4} = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global\AFR';
root{5} = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global\ASIA';
root{6} = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global\AUS';
root{7} = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global\CAM';
root{8} = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global\data\ARCTIC';

root_data = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global\data\data_prep\Population';

%% Load world population maps
disp('Loading Pop')
fname = sprintf('%s\\WordPop_ORNL_landscan2013.nc', root_data);
Pop_data = ncread(fname,'lspop2013')';
Seaid=Pop_data(1);
Pop_data(Pop_data==Seaid)=0;

%% Define continent area 4-corners
mult=60; %Multiplier to convert to 30s ((0.5*60)*60)/30=60)

% NAM
r1(1)=(61-1)*mult+1; %60 deg
r2(1)=130*mult; %25 deg
c1(1)=(85-1)*mult+1; %-138 deg
c2(1)=256*mult; %-52 deg
% EUR
r1(2)=(57-1)*mult+1; %62 deg
r2(2)=156*mult; %12 deg
c1(2)=(333-1)*mult+1; %-14 deg
c2(2)=500*mult; %70 deg
% SAM
r1(3)=(151-1)*mult+1; %15 deg
r2(3)=292*mult; %-56 deg
c1(3)=(175-1)*mult+1; %-93 deg
c2(3)=296*mult; %-32 deg
% AFR
r1(4)=(105-1)*mult+1; %38 deg
r2(4)=250*mult; %-35 deg
c1(4)=(323-1)*mult+1; %-19 deg
c2(4)=470*mult; %55 deg
% ASIA
r1(5)=(59-1)*mult+1; %61 deg
r2(5)=204*mult; %-12 deg
c1(5)=(475-1)*mult+1; %57 deg
c2(5)=720*mult; %180 deg
% AUS
r1(6)=(201-1)*mult+1; %-10 deg
r2(6)=292*mult; %-56 deg
c1(6)=(585-1)*mult+1; %112 deg
c2(6)=720*mult; %180 deg
% CAM
r1(7)=(103-1)*mult+1; %39 deg
r2(7)=170*mult; %5 deg
c1(7)=(123-1)*mult+1; %-119 deg
c2(7)=240*mult; %-60 deg
% ARCTIC
r1(8)=(1-1)*mult+1; %90 deg
r2(8)=60*mult; %60 deg
c1(8)=(1-1)*mult+1; %-180 deg
c2(8)=720*mult; %180 deg

%% Subset and resize data for continent
for rl=8
%     rl=7;
    fprintf('Continent #%d\n',rl)
    clear Pop_3 dir Pop nr10 nc10 nr nc
    
    disp('Selecting continent')
    Pop_3 = Pop_data(r1(rl):r2(rl),c1(rl):c2(rl));
    
    % figure(1);clf;imagesc(log(Pop_3));axis image;colormap(jet);
    
    %% Check alignment
    disp('Loading dir');
    fname = sprintf('%s\\dir.mat', root{rl});
    load(fname);
    [nr nc]=size(dir);
    
    %% Upscaler
    disp('Resizing')
    Pop = kron(Pop_3,int32(ones(1)));
    %Pop = Pop/4;
    
    %% Check
    [nr10 nc10]= size(Pop);

    assert(nr10==nr,'Size nr incorrect');
    assert(nc10==nc,'Size nc incorrect');
    
    %%
    disp('Saving')
    matfile = fullfile(root{rl}, sprintf('Pop.mat'));
    save(matfile,'Pop','-v7.3');
    
end


