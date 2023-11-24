clear all
close all
load('D:\GitConnect\data\ASIA\Basin_UIB\PantpeBasin_1.mat', 'acc')
load('D:\GitConnect\data\ASIA\Basin_UIB\PantpeBasin_1.mat', 'Q')
load('D:\GitConnect\data\ASIA\Basin_UIB\PantpeBasin_1.mat', 'outside')

acc(outside)=nan;
ma = nanmean(acc(:));
logma = 10^nanmean(log(acc(:)));

%% Plot acc
figure; 
subplot(1,3,1)
plot(sort(acc(:)))
set(gca, 'YScale', 'log')
yline(ma);
yline(logma);
title('acc sorted')
grid on
xlabel("cell ID")

subplot(1,3,2)
cdfplot(acc(:))
set(gca, 'XScale', 'log')
title("CDF of acc")

subplot(1,3,3)
plot(sort(acc(Q>=0.1)))
set(gca, 'YScale', 'log')
title("acc vals for Q>=0.1 cells")
xlabel("cell ID")
grid on

%% Get channels
channel1 = acc>ma & ~outside;
channel2 = acc>nanmedian(acc(:)) & ~outside;
channel3 = acc>(1/2*nanmean(acc(:))) & ~outside;
channel4 = acc>100 & ~outside;

%% of cells in diff channel defs
% sum(channel1(:))/sum(~isnan(acc(:)))*100
% sum(channel2(:))/sum(~isnan(acc(:)))*100
% sum(channel3(:))/sum(~isnan(acc(:)))*100
% sum(acc(:)>10)/sum(~isnan(acc(:)))*100
% sum(acc(:)>100)/sum(~isnan(acc(:)))*100
% sum(Q>=0.1 &~outside,'all' )/sum(~isnan(acc(:)))*100
% sum(Q>=0.1 & acc>logma &~outside,'all' )/sum(~isnan(acc(:)))*100

%% plots
figure;n=1;
subplot(3,3,n);n=n+1;
imagescnan(channel1)
title(sprintf("acc>mean(acc) or 866, %0.2f%% of basin cells",sum(channel1(:))/sum(~isnan(acc(:)))*100))

subplot(3,3,n);n=n+1;
imagescnan(channel3)
title(sprintf("acc>1/2*mean(acc), %0.2f%% of basin cells",sum(channel3(:))/sum(~isnan(acc(:)))*100))

subplot(3,3,n);n=n+1;
imagescnan(channel4)
title(sprintf("acc>100, %0.2f%% of basin cells",sum(channel3(:))/sum(~isnan(acc(:)))*100))

subplot(3,3,n);n=n+1;
imagescnan(Q>=0.1 & acc>logma &~outside)
title(sprintf("Q>0.1 & acc>logma, %0.2f%% of basin cells",sum(Q>=0.1 & acc>logma &~outside,'all' )/sum(~isnan(acc(:)))*100))

subplot(3,3,n);n=n+1;
imagescnan(Q>=0.1)
title(sprintf("Q>=0.1m^3/s, %0.2f%% of basin cells",sum(Q(:)>=0.1)/sum(~isnan(acc(:)))*100))

subplot(3,3,n);n=n+1;
imagescnan(acc>logma& ~outside)
title(sprintf("Acc>10^{mean(log(acc))} or 18.5, %0.2f%% of basin cells",sum(acc>logma& ~outside,'all')/sum(~isnan(acc(:)))*100))

subplot(3,3,n);n=n+1;
imagescnan(acc>10& ~outside)
title(sprintf("acc>10, %0.2f%% of basin cells",sum(acc(:)>10)/sum(~isnan(acc(:)))*100))

subplot(3,3,n);n=n+1;
imagescnan(channel2)
title(sprintf("acc>median(acc)), %0.2f%% of basin cells",sum(channel2(:))/sum(~isnan(acc(:)))*100))


%%
sum(channel3(:))/sum(~isnan(acc(:)))*100
sum(acc(:)>10)/sum(~isnan(acc(:)))*100
ma




