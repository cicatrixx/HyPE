% Create two figures showing change in total potential and distribution of
% potential with river segmentation
% Some figs from AGU 2020
clear all
close all
res="500m"; cellsz_m=500;

minQ=0.1; %m3/s
savefigs=0;
%subplot = @(m,n,p) subtightplot (m, n, p, [0.003 0.02]); %[vgap hgap]

%% Load dataset
root="G:/SurfDrive/GitConnect/output/fig_theoretical_fixedDEM";
load(sprintf("%s/%s/TheoreticalPot.mat",root,res));
cmap =linspecer(3);

datapath='G:/SurfDrive/GitConnect/data/UI/data';
data='UI500m_ArcGIS.mat';
fname=fullfile(datapath,data);
load(fname,'catchments','basinlabels', 'channel','outside')
%% Get sorted version of channel potentials and set <=0 to nan
dorangekm=dorange*cellsz_m/1000; % specified in terms of km
dorangename=strcat(string(dorangekm'),' km');

ncell=numel(cell_basinEnergy);
for i=1:length(dorangekm)
    sorted_Pot(:,i)=sort(reshape(channel_basinEnergy(:,:,i),[ncell,1]));
    sorted_Pot(sorted_Pot(:,i)<=0,i)=nan;
end

%% Get reshaped version
[nr, nc, nz]=size(channel_basinEnergy);
channel_re=reshape(channel_basinEnergy,[nr*nc,nz]);
channel_re(channel_re<=0)=nan;

%% Spatial plot of theory pot - very basic
figure
%imagescnan(~outside); hold all
imagescnan(maskBasin(channel_basinEnergy(:,:,1),~outside))
set(gca, 'ColorScale', 'log')
title("Theoretical potential at do=1(GWh)")
%savefig(sprintf("%s/%s/TheoreticalPotential_spatialatdo1.fig",root,res))


%% Plot matrix as scatter for do=1
sel_basinEnergy=channel_basinEnergy(:,:,1); %sel the do=1 case
idxo=  find(~isnan(sel_basinEnergy));
potval=sel_basinEnergy(idxo);
[ro,co]= ind2sub(size(sel_basinEnergy),idxo  );
[rch,cch]=find(channel);
nvalid=numel(ro);
inbasin=single(~outside);
inbasin(outside)=nan;

%
figure;
%imagescnan(inbasin);axis image;%colormap(gray)
imagescnan(maskBasin(catchments,~outside));axis image;%colormap(gray)
selbasins=1:8;
cmap=cbrewer('qual', 'Set3', length(basinlabels.basinIDs(selbasins)));
%cmap=linspecer(length(basinlabels.basinIDs(selbasins)),'qualitative');
colormap(cmap)
colorbar('Ticks', basinlabels.basinIDs(selbasins),'TickLabels',basinlabels.basinnames(selbasins),'Direction','reverse')

hold on
scatter(co,ro,rescale(potval,1,1500),'filled','MarkerFaceAlpha',0.4)%,'MarkerEdgeColor','blue')
%plot(cch,rch,'.b','markersize',.05,'DisplayName',"Channel") %plot channel
legend(sprintf("HP: %0.3g- %0.0f GWh",min(potval),max(potval)), "Channel")
title(sprintf("Channel Potenial: %0.2f TWh at %d sites",sum(potval,'omitnan')/1000,nvalid))
%h=colorbar('Location','eastoutside');
%h.Label.String='MW';

%% FINAL: Plot scatter points for range of do, TWh and lit vals
myval_allbasincells=sum(cell_basinEnergy,'all','omitnan'); %GWh
litvals_Indus=[channel_GWh(1) 	296368.32+523812.96	87600 201493	1944876]/1000; %in GWh converted to TWh
litvals_UI=[myval_allbasincells	296368.32+523812.96	87600  182908 	 1679935 ]/1000; %in GWh converted to TWh

litsrc={'Theoretical in current study'                                      ;
    'Visualized for Indus in India (CEA 2020) and Pakistan (WAPDA 2012)';
    'Theoretical for Pakistan (Lieftnick et al 1967)'                   ;
    'Technical for Upper Indus at 25km (Gernaat et al 2017)'              ;
    'Theoretical for Upper Indus at 0.22km (Hoes et al 2017)'               };

ccol = linspecer(length(litvals_UI),'qualitative');
figure
%subplot(1,2,1)
plot(dorangekm, channel_GWh/1000,"*-","LineWidth",2,'Color',.3*[1 1 1]);
hold all
yline(litvals_UI(1),"--","LineWidth",1.3,'Color',ccol(1,:)); %Visualized
yline(litvals_UI(2),":","LineWidth",1.3,'Color',ccol(2,:)); %Visualized
yline(litvals_UI(3),":","LineWidth",2,'Color',ccol(3,:)); %Lieftnick
yline(litvals_UI(4),":","LineWidth",2,'Color',ccol(4,:)); % Gernaat
yline(litvals_UI(5),":","LineWidth",1.3,'Color',ccol(5,:)); %Hoes

% Add text labels to ylines and one label for potential at min spacing
text(repmat(30,5,1),litvals_UI(1:5)*1.1,litsrc(1:5),'FontAngle','italic')
text(dorangekm(1),channel_GWh(1)/1000,sprintf('\\leftarrow %0.0f TWh/yr @ 500m', channel_GWh(1)/1000),'FontWeight','bold','Color',.3*[1 1 1]) % label max potential IDd
% add 4km label
text(dorangekm(5),channel_GWh(5)/1000,sprintf('\\leftarrow %0.0f TWh/yr @ 4km', channel_GWh(5)/1000),'FontWeight','bold','Color',.3*[1 1 1]) % label max potential IDd
% add 2km label
text(dorangekm(3),channel_GWh(3)/1000,sprintf('\\leftarrow %0.0f TWh/yr @ 2km', channel_GWh(3)/1000),'FontWeight','bold','Color',.3*[1 1 1]) % label max potential IDd

legend(litsrc(1))
%set(gca, 'YScale', 'log')
ylabel('Total energy (TWh/yr)')
xlabel('River segment length (km)')
title(sprintf("Total theoretical potential at %s and Q>=%0.1fm^3/s",res,minQ))
grid on

%% FINAL: Create boxplots for HP size distn varying w do
%figure;
subplot(1,2,2)
boxplot(channel_re,'Symbol','.')%,'Notch','on')

% add lines and labels for different sizes of HP -- from Siddiqui
HPclass=["Mega (>1000 MW)", "Large (500-1000 MW)", "Medium (50-500 MW)", "Small (5-50 MW)", "Mini (0.15-5 MW)", "Micro (0.005-0.15 MW)", "Pico (<0.005 MW)"];
HPsz_greaterthan=[1000	500	50	10	5	0.005]/1000*365*24; %MW converted to GWh assuming year round production

% from Hoes
%HPclass=["Large (>50MW)", "Small (5-50 MW)", "Mini (0.15-5 MW)", "Micro (0.005-0.15 MW)", "Pico (<0.005 MW)"];
%HPsz_greaterthan=[50, 5, 0.15, 0.005]*1000*365*24; %MW converted to kWh
hold on
for i=1:length(HPsz_greaterthan)
    yline(HPsz_greaterthan(i),'-.k');
end

text(repmat(size(channel_re,2),length(HPclass),1),[HPsz_greaterthan*5 HPsz_greaterthan(end)*.9] ,HPclass)
xticklabels(strcat(string(dorangekm'),'km'))
xlabel("River segment length")
ylabel("Theoretical potential at each segment (GWh)")
title(sprintf('Distribution of Energy in Valid Basin Cells at %s',res))
set(gca, 'YScale', 'log')
grid on
% savefig(fullfile(sprintf("%s/%s/TheoreticalPotential_channel.fig",root,res)))

%% Save fig as PDF
if savefigs
    savefig(sprintf("%s/%s/TheoreticalPotential_channel.fig",root,res))
    %orient('landscape')
    outpdf='G:\SurfDrive\GitConnect\output\fig_theoretical_fixedDEM\500m\PaperFig_total_distribution_wspacing.pdf';
    export_fig(outpdf)
    export_fig('G:\SurfDrive\GitConnect\output\fig_theoretical_fixedDEM\500m\PaperFig_total_distribution_wspacing.jpg')
end

%% FINAL: For 500m res Get subbasin-wise potential from rsi=1 channel potential
ch_basinEnergy=channel_basinEnergy(:,:,1);
selbasins=1:8;
cmap=linspecer(length(basinlabels.basinIDs(selbasins)),'qualitative');
for st=basinlabels.basinIDs(selbasins)'
    subPot_GWh(st-100) = sum(ch_basinEnergy(ch_basinEnergy > 0 & catchments ==st),'all');
    subPot_GWh_perarea(st-100) = subPot_GWh(st-100)/(sum(catchments ==st,'all')*cellsz_m^2);
end
subplot = @(m,n,p) subtightplot (m, n, p, [0.003 0.02]);

figure
%total potential
subplot(1,2,1)
%bar(subPot_GWh/1000)
bar(diag(subPot_GWh/1000),'stacked') % for colorful bars
ylabel('Total energy (TWh/yr)')
xticklabels(basinlabels.basinnames(selbasins))
xtickangle(45)
applymyplotformat('',cmap)

% per unit area
subplot(1,2,2)
%bar(subPot_GWh/1000)
bar(diag(subPot_GWh_perarea),'stacked') % for colorful bars
ylabel('Total energy per unit area (GWh/yr per m^2)')
xticklabels(basinlabels.basinnames(selbasins))
xtickangle(45)
applymyplotformat('',cmap)
% savefig(fullfile(sprintf("%s/%s/TheoreticalPotential_SubbasinWise.fig",root,res)))
% export_fig(fullfile(sprintf("%s/%s/TheoreticalPotential_SubbasinWise.jpg",root,res)))
figure
%subplot(2,2,3:4)
imagescnan(maskBasin(catchments,catchments>0&catchments<109))
colormap(cmap)
colorbar('Ticks', basinlabels.basinIDs(selbasins),'TickLabels',basinlabels.basinnames(selbasins),'Direction','reverse')
% savefig(fullfile(sprintf("%s/%s/TheoreticalPotential_SubbasinDefintion.fig",root,res)))

%% Tabulate total pot and number of sites in each potential class - all 6 classes
tblsum=table(flip(HPclass)','VariableNames',{'HPtype'});
tblcnt=tblsum;
for i = 1:length(dorangename)
    HPtype = discretize(channel_re(:,i),[0 flip(HPsz_greaterthan) Inf],'categorical',flip(HPclass));
    [potsum, cnt] = grpstats(channel_re(:,i),HPtype,{'sum','numel'});
    tblsum = addvars(tblsum,potsum,'NewVariableNames',dorangename(i));
    tblcnt = addvars(tblcnt,cnt,'NewVariableNames',dorangename(i));
    %tbltmp=grpstats(table(channel_re(:,i), HPtype), 'HPtype',{'sum'},'VarNames',{});
end

disp('% of small project (<50MW) in total')
prct_totpot=sum(tblsum{1:4,2:end})./sum(tblsum{:,2:end})*100
prct_count=sum(tblcnt{1:4,2:end})./sum(tblcnt{:,2:end})*100

%% FINAL: Barchart for distribution of potential - all classes
figure;
subplot(1,2,1)
bar(tblsum{:,2:end}','stacked')
%legend(tblsum.HPtype)
xticks(1:numel(dorangekm))
xticklabels(dorangename)
grid on
title("Total potential in each size category")
ylabel('Total energy (GWh/yr)')

subplot(1,2,2)
bar(tblcnt{:,2:end}','stacked')
legend(tblsum.HPtype)
xticks(1:numel(dorangekm))
xticklabels(dorangename)
ylabel('Number of sites')
grid on
title("Number of sites in each size category")
%set(gca, 'YScale', 'log')
% savefig(fullfile(sprintf("%s/%s/TheoreticalPotential_SizeWise.fig",root,res)))

%% GOOD: Proportion of small (<50MW) projects
figure
plot(prct_totpot,"*:")
hold all
plot(prct_count,".:")
grid on
legend('In terms of total potential','In terms of number of sites')
title('% of small project (<50MW) in total')
xticks(1:numel(dorangekm))
xticklabels(dorangename)
xtickangle(45)

ylabel('% ot total')
% savefig(fullfile(sprintf("%s/%s/TheoreticalPotential_PrctOfSmall.fig",root,res)))

%% GOOD: Barchart for distribution of potential - only subset distances
seldorange=[1 6 7 9 11 13 ];
figure;
subplot(1,2,1)
bar(tblsum{:,(seldorange+1)}')
legend(tblsum.HPtype)
xticks(1:numel(dorangekm(seldorange)))
xticklabels(dorangename(seldorange))
grid on
title("Total potential in each size category")
ylabel('Total potential energy at each site (GWh/yr)')

subplot(1,2,2)
bar(tblcnt{:,(seldorange+1)}')
xticks(1:numel(dorangekm(seldorange)))
xticklabels(dorangename(seldorange))
ylabel('Number of sites')
grid on
title("Number of sites in each size category")
set(gca, 'YScale', 'log')

%% Tabulate % and number of cells in each potential class - 2 classes
class2={'Small (<50MW)','Large'}';
tblsum2=table(class2,'VariableNames',{'HPtype'});
tblcnt2=tblsum2;
for i = 1:length(dorangename)
    HPtype = discretize(channel_re(:,i),[0 50/1000*365*24 Inf],'categorical',class2);
    [potsum, cnt] = grpstats(channel_re(:,i),HPtype,{'sum','numel'});
    tblsum2 = addvars(tblsum2,potsum,'NewVariableNames',dorangename(i));
    tblcnt2 = addvars(tblcnt2,cnt,'NewVariableNames',dorangename(i));
    %tbltmp=grpstats(table(channel_re(:,i), HPtype), 'HPtype',{'sum'},'VarNames',{});
end