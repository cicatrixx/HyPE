% Read outputs of EvalTheorticalPot for the three resolutions and compile
% them - only at Upper Indus level
% Last Updated: 22 Sept for David meeting

clear all
close all
res=["15sUI", "500m", "5km"];%"15s"
Qplots=0;       %Prepare Q comparison figs
showplot=1;     %Show tmp fig

%% Load all three datasets
root="G:/SurfDrive/GitConnect/output/figs_theoretical/";
b15s=load(sprintf("%s/%s/TheoreticalPot.mat",root,res{1}));
b5km=load(sprintf("%s/%s/TheoreticalPot.mat",root,res{3}));
b500m_old=load(sprintf("%s/%s/TheoreticalPot.mat",root,res{2}));

b500m=load(sprintf("G:/SurfDrive/GitConnect/output/fig_theoretical_fixedDEM/%s/TheoreticalPot_olddo.mat",res{2}));
cmap =linspecer(3);

%% FINAL: Compare pot vs do for different datasets Plot range of do in km, potential in TWh + litvals
%dorange=[1,3,5,10:10:100];
dorangekm=ceil([500, 1000, 2000, 3000, 4000, 5000, 10000:10000:100000])/1000; % specified in terms of km
figure
plot(dorangekm,b15s.channel_GWh/1000,".:","Color",cmap(1,:))
hold all
plot(dorangekm(6:end),b5km.channel_GWh(6:end)/1000,".:","Color",cmap(3,:)) % skip dorange<5km
plot(dorangekm,b500m_old.channel_GWh/1000,".:","Color",cmap(2,:))
plot(dorangekm,b500m.channel_GWh/1000,".-","Color",cmap(2,:))

legend(["15sUI", "5km", "500m-old", "500m-new"])
ylabel('Total energy (TWh/yr)')
xlabel('River spacing (km)')
applymyplotformat("Potential at 3 resolutions")
%set(gca, 'YScale', 'log')

%% add litvals
litvals_UI=[0 296368.32+523812.96	87600  182908 	 1679935 ]/1000; %in GWh converted to TWh
litsrc={      'tmp'                              ;
    'Visualized for Indus in India (CEA 2020) and Pakistan (WAPDA 2012)';
    'Theoretical for Pakistan (Lieftnick et al 1967)'                   ;
    'Technical for Upper Indus at 25km (Gernaat et al 2017)'              ;
    'Theoretical for Upper Indus at 0.22km (Hoes et al 2017)'               };

ccol =[1 1 1]*0.2; %linspecer(length(litvals_UI),'qualitative');
yline(litvals_UI(2),"--","LineWidth",1.2,'Color',ccol); %Visualized
yline(litvals_UI(3),"-","LineWidth",1.1,'Color',ccol); %Lieftnick
yline(litvals_UI(4),":","LineWidth",2,'Color',ccol); % Gernaat
yline(litvals_UI(5),"-","LineWidth",1.1,'Color',ccol); %Hoes
text(repmat(100,4,1),litvals_UI(2:5)*1.1,litsrc(2:5))

%% FINAL: Distribution of Energy in Valid Basin Cells - invalid cells are already nan
figure
plot(sort(b15s.cell_basinEnergy(b15s.cell_basinEnergy>0)),"Color",cmap(1,:)) % only plot non-zeros
hold all
plot(sort(b500m.cell_basinEnergy(b500m.cell_basinEnergy>0)),"Color",cmap(2,:)) % only plot non-zeros
plot(sort(b5km.cell_basinEnergy(b5km.cell_basinEnergy>0)),"Color",cmap(3,:)) % only plot non-zeros
grid on
axis square
set(gca, 'YScale', 'log', 'XScale', 'log')
xlabel("Cell number sorted by cell potential")
ylabel("Theoretical Potential at Cell (GWh)")
title('Distribution of Energy in Valid Basin Cells')
legend(res)

%% Output: Plot range of do in km, potential in GWh + theoretical potentials
if showplot
    %% Plot range of do in km, potential in GWh + theoretical potential
    figure
    yline(sum(b15s.cell_basinEnergy(b15s.cell_basinEnergy(:)>=0)) ,"--","Color",cmap(1,:));
    hold all
    plot(dorangekm,b15s.channel_GWh,".-","Color",cmap(1,:))
    yline(sum(b500m.cell_basinEnergy(b500m.cell_basinEnergy(:)>=0)) ,"--","Color",cmap(2,:));
    plot(dorangekm,b500m.channel_GWh,".-","Color",cmap(2,:))
    yline(sum(b5km.cell_basinEnergy(b5km.cell_basinEnergy(:)>=0)) ,"--","Color",cmap(3,:));
    plot(dorangekm(6:end),b5km.channel_GWh(6:end),".-","Color",cmap(3,:)) % skip dorange<5km
    
    %
    ylabel('Total channel energy (GWh)')
    xlabel('River spacing (km)')
    title(sprintf("Cell-by-cell vs Channel Potential at 3 resolutions"))
    grid on
    lgd=legend([repmat(" ",4,1);'Cell by cell potential'; 'Channel potential'], 'Location', 'Best','NumColumns', 3);
    lgd.Title.String=sprintf("%s                 %s         %s                          Type",res);
    box on
    %set(gca, 'YScale', 'log')
    
    %% Output: Plot range of do in km, potential as % of theoretical
    figure
    plot(dorangekm,b15s.channel_GWh*100/sum(b15s.cell_basinEnergy(b15s.cell_basinEnergy(:)>=0)),".-","Color",cmap(1,:))
    hold all
    plot(dorangekm,b500m.channel_GWh*100/sum(b500m.cell_basinEnergy(b500m.cell_basinEnergy(:)>=0)),".-","Color",cmap(2,:))
    plot(dorangekm(6:end),b5km.channel_GWh(6:end)*100/sum(b5km.cell_basinEnergy(b5km.cell_basinEnergy(:)>=0)),".-","Color",cmap(3,:))
    ylabel('% of cell by cell theoretical energy')
    xlabel('River spacing (km)')
    title(sprintf("Cell-by-cell vs Channel Potential at 3 resolutions"))
    grid on
    legend(res)
    box on
    
    %set(gca, 'YScale', 'log')
    
    %% Output: Plot range of do, GWh as % of theoretical in cell spacing
    figure
    plot(dorangekm/.450,b15s.channel_GWh*100/sum(b15s.cell_basinEnergy(b15s.cell_basinEnergy(:)>=0)),".-","Color",cmap(1,:))
    hold all
    plot(dorangekm/.500,b500m.channel_GWh*100/sum(b500m.cell_basinEnergy(b500m.cell_basinEnergy(:)>=0)),".-","Color",cmap(2,:))
    plot(dorangekm/5,b5km.channel_GWh*100/sum(b5km.cell_basinEnergy(b5km.cell_basinEnergy(:)>=0)),".-","Color",cmap(3,:))
    ylabel('% of theoretical energy')
    xlabel('Channel spacing in # of cells')
    title(sprintf("Channel Potential as percentage of cell-by-cell potential at 3 resolutions"))
    grid on
    legend(res)
    
    %% cdf of Energy in Valid Basin Cells
    figure
    cdfplot(b15s.cell_basinEnergy(b15s.cell_basinEnergy>0)) %,"Color",cmap(1,:)) % only plot non-zeros
    hold all
    cdfplot(b500m.cell_basinEnergy(b500m.cell_basinEnergy>0))%,"Color",cmap(2,:)) % only plot non-zeros
    cdfplot(b5km.cell_basinEnergy(b5km.cell_basinEnergy>0))%,"Color",cmap(3,:)) % only plot non-zeros
    grid on
    axis square
    set(gca, 'XScale', 'log')%, 'YScale', 'log'
    %ylabel("Theoretical Potential at Cell (GWh)")
    title('CDF of Energy in Valid Basin Cells')
    legend(res)
    colororder(cmap(1:3,:))
end

%% Q comparison plots
if Qplots
    %% Load Qs
    fname='G:/PaperData/DavidGernaat/Sanita_model_package/data/ASIA/Basin/Basin_5.mat';
    load(fname,'Q');
    b15s.Q=Q;
    b15s.Q(isnan(b15s.cell_Hgross(:)))=nan;
    
    %
    QMLAB="MLAB500m_40yrClimatology_m3day.mat";
    load(fullfile(datapath,QMLAB),'Q500m_m3day');
    b500m.Q=Q500m_m3day{13}/(24*60*60);
    b500m.Q(isnan(b500m.cell_Hgross(:)))=nan;
    %
    QSPHY = 'SPHY5km_40yrClimatology.mat';
    load(fullfile(datapath,QSPHY),'Qm_SPHY_m3sec');
    b5km.Q=Qm_SPHY_m3sec{13};
    b5km.Q(isnan(b5km.cell_Hgross(:)))=nan;
    %
    
    %% Plot overlaid cell Hgross and Q for each resolutions in one subfigure
    figure
    i=1;
    [~,b15s.order]=sort(b15s.cell_Hgross(:));
    subplot(1,3,i)
    plot(b15s.cell_Hgross(b15s.order),"Color",cmap(i,:),"LineWidth",1.5);
    hold all
    plot(b15s.Q(b15s.order),":","Color",cmap(i,:))
    grid on
    legend("Hgross in m", "Discharge in m^3/s","Location","Best")
    title(res(i))
    %
    i=2;
    [~,b500m.order]=sort(b500m.cell_Hgross(:));
    subplot(1,3,i)
    plot(b500m.cell_Hgross(b500m.order),"Color",cmap(i,:),"LineWidth",1.5);
    hold all
    plot(b500m.Q(b500m.order),":","Color",cmap(i,:))
    grid on
    legend("Hgross in m", "Discharge in m^3/s","Location","Best")
    title(res(i))
    
    %
    i=3;
    [~,b5km.order]=sort(b5km.cell_Hgross(:));
    subplot(1,3,i)
    plot(b5km.cell_Hgross(b5km.order),"Color",cmap(i,:),"LineWidth",1.5);
    hold all
    plot(b5km.Q(b5km.order),":","Color",cmap(i,:))
    grid on
    legend("Hgross in m", "Discharge in m^3/s","Location","Best")
    title(res(i))
    
    %% Plot overlaid cell Energy, Hgross and Q for each resolutions in one subfigure
    %
    figure
    i=1;
    [~,b15s.order]=sort(b15s.cell_basinEnergy(:));
    subplot(1,3,i)
    plot(b15s.cell_Hgross(b15s.order),"Color",cmap(i,:));
    hold all
    plot(b15s.Q(b15s.order),":k")%,"Color",cmap(i,:))
    plot(b15s.cell_basinEnergy(b15s.order)/10)%,"Color",cmap(i,:))
    grid on
    legend("Hgross in m", "Discharge in m^3/s","Location","Best")
    title(res(i))
    %
    i=2;
    [~,b500m.order]=sort(b500m.cell_basinEnergy(:));
    subplot(1,3,i)
    plot(b500m.cell_Hgross(b500m.order),"Color",cmap(i,:));
    hold all
    plot(b500m.Q(b500m.order),":k")%,"Color",cmap(i,:))
    plot(b500m.cell_basinEnergy(b500m.order)/10)%,"Color",cmap(i,:))
    grid on
    legend("Hgross in m", "Discharge in m^3/s","Location","Best")
    title(res(i))
    
    %
    i=3;
    [~,b5km.order]=sort(b5km.cell_basinEnergy(:));
    subplot(1,3,i)
    plot(b5km.cell_Hgross(b5km.order),"Color",cmap(i,:));
    hold all
    plot(b5km.Q(b5km.order),":k")%"Color",cmap(i,:))
    plot(b5km.cell_basinEnergy(b5km.order)/10)%,"Color",cmap(i,:))
    grid on
    legend("Hgross in m", "Discharge in m^3/s","Location","Best")
    title(res(i))
end