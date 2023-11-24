% figs for SustaIndus Project Group meeting on 30 Sept 2020
% Eval # of prrojects and size under each HP size category
% Compile binary layers
% Eval reduction in potential w binary count

clear all
close all
res=["15sUI","500m","5km"];%"15s"

showplot=1;  %Show tmp fig
dorangekm=ceil([500, 1000, 2000, 3000, 4000, 5000, 10000:10000:100000])/1000; % specified in terms of km
dorangename=strcat(string(dorangekm'),' km');

%add labels for different sizes of HP
HPclass=["Mega (>1000 MW)", "Large (500-1000 MW)", "Medium (50-500 MW)", "Small-II (10-50 MW)", "Small-I (5-10 MW)", "Mini (0.15-5 MW)", "Micro (<0.15 MW)"];
HPsz_greaterthan=[1000	500	50	10	5	0.15]/1000*365*24; %MW converted to GWh

%% Load all three datasets
root="G:/SurfDrive/GitConnect/output/figs_theoretical/";
b15s=load(sprintf("%s/%s/TheoreticalPot.mat",root,res{1}));
b500m=load(sprintf("%s/%s/TheoreticalPot.mat",root,res{2}));
b5km=load(sprintf("%s/%s/TheoreticalPot.mat",root,res{3}));
cmap =linspecer(3);

%% Get sorted version of basin cell potentials and set <=0 to nan
mydata=b500m;
ncell=numel(mydata.cell_basinEnergy);
for i=1:length(dorangekm)
    sorted_Pot(:,i)=sort(reshape(mydata.channel_basinEnergy(:,:,i),[ncell,1]));
    sorted_Pot(sorted_Pot(:,i)<=0,i)=nan;
end

%% Get reshaped version
[nr, nc, nz]=size(mydata.channel_basinEnergy);
channel_re=reshape(mydata.channel_basinEnergy,[nr*nc,nz]);
channel_re(channel_re<=0)=nan;

%% GOOD: Create boxplots for HP size distn varying w do
figure;boxplot(channel_re,'Symbol','.')%,'Notch','on')

hold on
yline(HPsz_greaterthan(1),'-.k');
yline(HPsz_greaterthan(2),'-.k');
yline(HPsz_greaterthan(3),'-.k');
yline(HPsz_greaterthan(4),'-.k');
yline(HPsz_greaterthan(5),'-.k');
yline(HPsz_greaterthan(6),'-.k');
text(repmat(17,7,1),[HPsz_greaterthan*5 HPsz_greaterthan(end)*.9] ,HPclass)

xticklabels(strcat(string(dorangekm'),'km'))
xlabel("River Spacing")
ylabel("Theoretical Potential at Cell (GWh)")
title(sprintf('Distribution of Energy in Valid Basin Cells at %s',res{2}))
set(gca, 'YScale', 'log')
grid on


%% Tabulate % and number of cells in each potential class - all 6 classes
tblsum=table(flip(HPclass)','VariableNames',{'HPtype'});
tblcnt=tblsum;
for i = 1:length(dorangename)
    HPtype = discretize(channel_re(:,i),[0 flip(HPsz_greaterthan) Inf],'categorical',flip(HPclass));
    [sum, cnt] = grpstats(channel_re(:,i),HPtype,{'sum','numel'});
    tblsum = addvars(tblsum,sum,'NewVariableNames',dorangename(i));        
    tblcnt = addvars(tblcnt,cnt,'NewVariableNames',dorangename(i));
    %tbltmp=grpstats(table(channel_re(:,i), HPtype), 'HPtype',{'sum'},'VarNames',{});
end
% plot for last dorange
figure;histogram(HPtype) 

%% Tabulate % and number of cells in each potential class - 2 classes
class2={'Small (<50MW)','Large'}';
tblsum2=table(class2,'VariableNames',{'HPtype'});
tblcnt2=tblsum2;
for i = 1:length(dorangename)
    HPtype = discretize(channel_re(:,i),[0 50/1000*365*24 Inf],'categorical',class2);
    [sum, cnt] = grpstats(channel_re(:,i),HPtype,{'sum','numel'});
    tblsum2 = addvars(tblsum2,sum,'NewVariableNames',dorangename(i));        
    tblcnt2 = addvars(tblcnt2,cnt,'NewVariableNames',dorangename(i));
    %tbltmp=grpstats(table(channel_re(:,i), HPtype), 'HPtype',{'sum'},'VarNames',{});
end

%% Create histogram plots - bad overlap
%edges(1) is the left edge of the first bin, and edges(end) is the right edge of the last bin.
%The value X(i) is in the kth bin if edges(k) ? X(i) < edges(k+1)
bins4HPsize_name={'Micro';'Mini';'Small1';'Small2';'Medium'; 'Large';	'Mega'};
bins4HPsize = [0.0001 1.314	43.8	87.6	438	4380	8760 ceil(max(sorted_Pot,[],'all'))]; % in GWh/yr
figure
i=1;
h(i)=histogram(sorted_Pot(:,i),bins4HPsize);
hold on
for i=2:length(dorangekm)
    %figure
    h(i)=histogram(sorted_Pot(:,i),bins4HPsize);
    set(gca, 'XScale', 'log')
end

%% GOOD: Try stackedplot stairs - not bad but not the best
% get counts in each HP size
for i=1:length(dorangekm)
    bincounts(:,i) = histcounts(sorted_Pot(:,i),bins4HPsize);
end

figure; s=stackedplot(1:length(bins4HPsize)-1,bincounts);
s.LineWidth=1;
for i=1:length(dorangekm)
    s.LineProperties(i).PlotType='stairs';
end

tbl=array2table(bincounts,'RowNames',bins4HPsize_name,'VariableNames',strcat('@',string(dorangekm'),' km'));
figure; s=stackedplot(tbl);
xticklabels(bins4HPsize_name)

%% Get cell counts in each binary layer
theoryP=b500m.channel_basinEnergy(:,:,1);
theoryP(theoryP<=0)=nan;
theory_ncell=sum(~isnan(theoryP),'all');
theoryPtot=b500m.channel_GWh(1);
figure;imagescnan(theoryP);axis image
set(gca, 'ColorScale', 'log')

%% Load layers
x(1)=load('G:\SurfDrive\GitConnect\data\UI\data\WaterBodies.mat');
x(2)=load('G:\SurfDrive\GitConnect\data\UI\data\Glaciers.mat');
x(3)=load('G:\SurfDrive\GitConnect\data\UI\data\MountainAreas.mat');
x(4)=load('G:\SurfDrive\GitConnect\data\UI\data\WDPA.mat');
x(5)=load('G:\SurfDrive\GitConnect\data\UI\data\GSHAP.mat');
x(6)=load('G:\SurfDrive\GitConnect\data\UI\data\LandSlideSusceptibility.mat');

%% high risk hazards
x(5).data=x(5).data>=4.0; % seismic>4pga
x(6).data=x(6).data>=4.0; % landslide

%%
for i=1:6
    figure;imagescnan(x(i).data>0)
    valid(:,:,i) = x(i).data>0;
    Exclude_ncell(i,1) = sum(valid(:,:,i) & ~isnan(theoryP),'all');
    Exclude_P(i,1) = sum(theoryP(valid(:,:,i) & ~isnan(theoryP)),'all');
end
Exclude_ncell_prct=Exclude_ncell/theory_ncell*100;
 Exclude_P_prct=Exclude_P/theoryPtot*100;

tbl=table(Exclude_ncell,Exclude_ncell_prct , Exclude_P,Exclude_P_prct)

%% Distance to transmission lines
load('G:\SurfDrive\GitConnect\data\UI\data\DisCost.mat')
celldist2line=data;
cellsz=0.5; %in km
celldist2line(celldist2line==0)=nan;

figure;imagescnan(celldist2line)
figure;boxplot(celldist2line(:)*cellsz,'Symbol','.')




