% Compare Q TS averaging to monthly LT avgs in historical
% This verifies that my monTS data is indeed in m3/month and needs to be
% changed to compare w the LTA in monLTA.

% There is -7 to 7% difference between the Q_LTA evaluated here from Q_monTS and the
% Q_monLTA evaluated from monLTA runoff files. That may be reasonable w
% rounding errors and errors due to rescaling or R to Q.

clc
clear all
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
root=pwd;
matpath=fullfile(root,"data","UI","data");
%% Get fpaths to Q_TS files
fpaths_TS=path2fldrfiles(matpath,"MLAB500m_FutQ_TS_m3day_*");
fpaths_monLT=[path2fldrfiles(matpath,"MLAB500m_40yr*")]; % path2fldrfiles(matpath,"MLAB500m_FutQ_m3day_*")];

%% Get basin bounds
disp('Loading acc and dir at 500m');
data500m='UI500m_ArcGIS.mat';
load(fullfile(matpath,data500m), 'outside', 'outlets','basinlabels')%,'basinmask','channel') %,  'fdir'

%% load Q_TS and QLTA
f=1;
load(fpaths_TS{f},'Q500m_m3day');
mymonTS_m3month=cat(3,Q500m_m3day{:}); 
Q_LTA_m3_day=load(fpaths_monLT{f},'Q500m_m3day');
Q_m3s_SPHY=maskBasin(cat(3,Q_LTA_m3_day.Q500m_m3day{:})/(24*60*60),~outside); % convert cell to matrix

%% Eval and plot Qdiff for hist TS and hist LTA
m=12;
cTS=maskBasin(mymonTS_m3month(:,:,m+12*22)/eomday(1,m),~outside);
cLTA=maskBasin(Q_LTA_m3_day.Q500m_m3day{m}, ~outside);
figure;
subplot(2,2,1)          
imagescnan(cTS)
colorbar
set(gca,'colorscale','log')
title(sprintf('Q TS for m=%d', m))

subplot(2,2,2)
imagescnan(cLTA)
colorbar
title(sprintf('Q LTA for m=%d', m))
set(gca,'colorscale','log')

subplot(2,2,3)
diffprct=(cTS-cLTA)./cLTA*100;
imagescnan(diffprct)
caxis([-10 100])
title(sprintf('Q TS-LTA in %% for m=%d', m))
colorbar

subplot(2,2,4)
boxchart(diffprct(:))%,'MarkerStyle','none')
title(sprintf('Q TS-LTA in %% for m=%d', m))

sgtitle(["For", fpaths_TS{f}], 'Interpreter','none')

%% Extract station TS from Q500m
%find all outlets within basin
outletID=unique(outlets(outlets>0)); %~isnan(outlet)
npts=length(outletID);
for pt=1:8
    % Get 500m TS at outlets
    [r0,c0]=find(outlets==outletID(pt));
    Q500m_TS_compiled(:,pt) = squeeze(mymonTS_m3month(r0,c0,:));
    Q500m_LTA_old(1:13,pt) = squeeze(Q_m3s_SPHY(r0,c0,:));
end

%% Create timeseries for days in month
hist_yrs=40;
hist_ts = datetime('January 1, 1979') + calmonths(0:12*hist_yrs-1)';
hist_daysinmon_ts=eomday(year(hist_ts),month(hist_ts));
myTS_m3s= Q500m_TS_compiled./hist_daysinmon_ts/(24*60*60);

%% Get LTA from TS
myLTA=plotMonthlyAverageDaily(myTS_m3s',0,1,basinlabels.basinnames,"Discharge (m^3/s)");

% Add old LTA
for pt=1:8
subplot(4,2,pt)
plot(Q500m_LTA_old(1:12,pt),'xk','LineWidth',8,'DisplayName',"Hydrus LTA")
end

%% Get diff
diffprct=(myLTA'-Q500m_LTA_old(1:12,:))./Q500m_LTA_old(1:12,:)*100;
figure
plot(diffprct)
title("% difference between myLTA_Q from monTS_Q vs LTA_Q from LTA runoffs")
[max(diffprct,[],'all') min(diffprct,[],'all')]
colororder(hsv(8))
legend(basinlabels.basinnames)
grid on