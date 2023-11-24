%% Compare cost vs size across the portfolios for different scenarios

% load all tech-tech and sust-tech portfolios then filter and plot only tech-fin and sust-fin
% basin and sub-basin wise mean and median plots

clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')
ofmatname=fullfile(rootof,'MainFutScenarios.mat');

load(ofmatname, 'pcsout',...
    'runnames','basindata','catchments_cl')

%% Compile all scenarios and their details
runnames_cl=strrep(strrep(runnames,'_Tech_Fin','_Technical'),'_Sust_RiskAverse','_Sustainable')';

fdetails=split(runnames_cl,'_');
fdetails(1:2,1)={'Historical';'Historical'};
fdetails(1:2,2)={'ERA5-Land';'ERA5-Land'};
fdetails(1:2,3)={'1985';'1985'};
fdetails(1:2,4)={'2014';'2014'};

tflist=join(fdetails(:,3:4),"-");
cornerlist=categorical(['NA';'NA'; repmat(repelem(cornernames,2),1,6)'],['NA', cornernames]);
scen_details = table(categorical(fdetails(:,1),unique(fdetails(:,1))),cornerlist,fdetails(:,2),categorical(tflist,unique(tflist)),categorical(fdetails(:,5),unique(fdetails(:,5))), ...
    'VariableNames',{'ssp','corner','gcm','yrs','pottype'});
%scen_details = table(categorical(fdetails(:,1),unique(fdetails(:,1)),'Ordinal',1),cornerlist,fdetails(:,2),categorical(tflist,unique(tflist),'Ordinal',1),fdetails(:,5), ...

%% Compile all scenario portfolios into one big table data
scen_data=table();
for scn=1:50
    scen_data_tmp=struct2table(rmfield(pcsout{scn},{'PID','co_arc','ro_arc','runname'}));
    %table(pcsout{scn}.PnetAlls,pcsout{scn}.COEAlls, pcsout{scn}.SysIDs,pcsout{scn}.subbasinIDs,...
    % 'VariableNames',{'GWh/yr','USD/kWh','PlantType','Subbasin'});

    %remove rows with missing data
    scen_data_tmp=rmmissing(scen_data_tmp);
    scen_data_tmp(:,'runname')={pcsout{scn}.runname};
    scen_data_tmp=[scen_data_tmp repmat(scen_details(scn,:),height(scen_data_tmp),1)];
    scen_data=[scen_data; scen_data_tmp];
end

scen_data.subname= categorical(scen_data.subbasinIDs,basindata.basinIDs,basindata.basinnames);
scen_data.planttype= categorical(scen_data.SysIDs,[2 1], planttypes);

disp("Created long table w HP data")

selfinpot=scen_data.COEAlls<=costlim;
%% Box plot for GWh - Basin level
figure
%subplot(2,1,1)
boxplot(scen_data{selfinpot,"PnetAlls"},scen_data{selfinpot,["pottype","planttype","yrs","ssp","corner"] },'whisker',inf,...
    'ColorGroup',scen_data{selfinpot,"corner"},'Colors',[0 0 0; cbrewer2("Dark2",4)],...%'orientation','horizontal',...
    'PlotStyle','compact','OutlierSize',2,'Symbol','+','FactorDirection','data','FullFactors','off','FactorSeparator',[1,2]); %'FactorGap',[10 ,.1],
grid on

% using inf whiskers or no only affects if outliers are seen and not the
% median values
%'whisker',inf,'color',c_ssp_main,
%set(gca,'YScale', 'log','XGrid','off','FontSize',10) %,'XTicklabel',[])
ylabel('Project Energy (GWh per year)','FontWeight','bold')
%xticklabels([])

%% Box plot for COE - Basin level
figure
%subplot(2,1,2)
boxplot(scen_data{selfinpot,"COEAlls"},scen_data{selfinpot,["pottype","planttype","yrs","ssp","corner"] },'whisker',inf,...
    'ColorGroup',scen_data{selfinpot,"corner"},'Colors',[0 0 0; cbrewer2("Dark2",4)],...%'orientation','horizontal',...
    'PlotStyle','compact','OutlierSize',2,'Symbol','+','FactorDirection','data','FullFactors','off','FactorSeparator',[1,2]); %'FactorGap',[10 ,.1],
grid on
%'whisker',inf,'color',c_ssp_main,
%set(gca,'YScale', 'log','XGrid','off','FontSize',10) %,'XTicklabel',[])
ylabel('Project Per unit production cost (USD2010 per KWh)','FontWeight','bold')

%% Box plot for GWh - SubBasin level
for p=pottypes2
    selrows=scen_data.pottype==p;
    figure
    %subplot(2,1,1)
    boxplot(scen_data{selrows,"PnetAlls"},scen_data{selrows,["subname","yrs","ssp","corner"] },'whisker',inf,'ColorGroup',scen_data{selrows,"corner"},'Colors',[0 0 0; cbrewer2("Dark2",4)],...%'orientation','horizontal',...
        'PlotStyle','compact','OutlierSize',2,'Symbol','+','FactorDirection','data','FullFactors','off','FactorSeparator',[1,2]); %'FactorGap',[10 ,.1],
    grid on
    title(p)
    % using inf whiskers or no only affects if outliers are seen and not the
    % median values
    %'whisker',inf,'color',c_ssp_main,
    set(gca,'YScale', 'log','XGrid','off','FontSize',10) %,'XTicklabel',[])
    ylabel('Project Energy (GWh per year)','FontWeight','bold')
    %xticklabels([])
end

%% Box plot for GWh - SubBasin level ONLY FINANCIAL
for p=pottypes2(1)
    selrows=scen_data.pottype==p & scen_data.COEAlls<=costlim;
    boxmean=groupsummary(scen_data{selrows,"PnetAlls"},scen_data{selrows,["subname","yrs","ssp","corner"]},'mean');
    figure
    %subplot(2,1,1)
    boxplot(scen_data{selrows,"PnetAlls"},scen_data{selrows,["subname","yrs","ssp","corner"] },'whisker',inf,'ColorGroup',scen_data{selrows,"corner"},'Colors',[0 0 0; cbrewer2("Dark2",4)],...%'orientation','horizontal',...
        'PlotStyle','compact','OutlierSize',2,'Symbol','+','FactorDirection','data','FullFactors','off','FactorSeparator',[1,2]); %'FactorGap',[10 ,.1],
    hold on
    plot(boxmean,'*')
    grid on
    title(strcat(p,'-Financial'))
    % using inf whiskers or no only affects if outliers are seen and not the
    % median values
    %'whisker',inf,'color',c_ssp_main,
    set(gca,'YScale', 'log','XGrid','off','FontSize',10) %,'XTicklabel',[])
    ylabel('Project Energy (GWh per year)','FontWeight','bold')
    %xticklabels([])
end

%% Box plot for COE - SubBasin level
for p=pottypes2
    selrows=scen_data.pottype==p;
    figure
    %subplot(2,1,1)
    boxplot(scen_data{selrows,"COEAlls"},scen_data{selrows,["subname","yrs","ssp","corner"] },'whisker',inf,'ColorGroup',scen_data{selrows,"corner"},'Colors',[0 0 0; cbrewer2("Dark2",4)],...%'orientation','horizontal',...
        'PlotStyle','compact','OutlierSize',2,'Symbol','+','FactorDirection','data','FullFactors','off','FactorSeparator',[1,2]); %'FactorGap',[10 ,.1],
    grid on
    title(p)
    % using inf whiskers or no only affects if outliers are seen and not the
    % median values
    %'whisker',inf,'color',c_ssp_main,
    set(gca,'YScale', 'log','XGrid','off','FontSize',10) %,'XTicklabel',[])
    ylabel('Project Per unit production cost (USD2010 per KWh)','FontWeight','bold')
    %xticklabels([])
end

%% Get medians and means for all or only financial project
onlyfinprjs=01;
if ~onlyfinprjs
    % this one would be technical pot
    coe_basin_median = groupsummary(scen_data,["pottype","planttype","yrs","ssp","corner"],"median","COEAlls");
    coe_basin_mean = groupsummary(scen_data,["pottype","planttype","yrs","ssp","corner"],"mean","COEAlls");

    gwh_basin_median = groupsummary(scen_data,["pottype","planttype","yrs","ssp","corner"],"median","PnetAlls");
    gwh_basin_mean = groupsummary(scen_data,["pottype","planttype","yrs","ssp","corner"],"mean","PnetAlls");
    %gwh_basin_mean = groupsummary(scen_data.PnetAlls,scen_data{:,["pottype","planttype","yrs","ssp","corner"]},"mean");


    coe_subbasin_median = groupsummary(scen_data,["subname","pottype","planttype","yrs","ssp","corner"],"median","COEAlls");
    coe_subbasin_mean = groupsummary(scen_data,["subname","pottype","planttype","yrs","ssp","corner"],"mean","COEAlls");

    gwh_subbasin_median = groupsummary(scen_data,["subname","pottype","planttype","yrs","ssp","corner"],"median","PnetAlls");
    gwh_subbasin_mean = groupsummary(scen_data,["subname","pottype","planttype","yrs","ssp","corner"],"mean","PnetAlls");
    ff=1; % idx for separating figs
else
    mydata=scen_data(selfinpot,:);
    %rename pot types to indicat it is only financial
    mydata.pottype=renamecats(mydata.pottype,pottypes2,strcat(pottypes2,"-Fin"));

    coe_basin_median = groupsummary(mydata,["pottype","planttype","yrs","ssp","corner"],"median","COEAlls");
    coe_basin_mean = groupsummary(mydata,["pottype","planttype","yrs","ssp","corner"],"mean","COEAlls");

    gwh_basin_median = groupsummary(mydata,["pottype","planttype","yrs","ssp","corner"],"median","PnetAlls");
    gwh_basin_mean = groupsummary(mydata,["pottype","planttype","yrs","ssp","corner"],"mean","PnetAlls");
    %gwh_basin_mean = groupsummary(scen_data.PnetAlls,scen_data{:,["pottype","planttype","yrs","ssp","corner"]},"mean");


    coe_subbasin_median = groupsummary(mydata,["subname","pottype","planttype","yrs","ssp","corner"],"median","COEAlls");
    coe_subbasin_mean = groupsummary(mydata,["subname","pottype","planttype","yrs","ssp","corner"],"mean","COEAlls");

    gwh_subbasin_median = groupsummary(mydata,["subname","pottype","planttype","yrs","ssp","corner"],"median","PnetAlls");
    gwh_subbasin_mean = groupsummary(mydata,["subname","pottype","planttype","yrs","ssp","corner"],"mean","PnetAlls");
    ff=10; % idx for separating figs
end

%% Get medians for technical potential only
selrows=scen_data.pottype=="Technical";
mydata=scen_data(selrows,:);
coe_basin_median_tech = groupsummary(mydata,["pottype","planttype","yrs","ssp","corner"],"median","COEAlls");
coe_basin_mean_tech = groupsummary(mydata,["pottype","planttype","yrs","ssp","corner"],"mean","COEAlls");

gwh_basin_median_tech = groupsummary(mydata,["pottype","planttype","yrs","ssp","corner"],"median","PnetAlls");
gwh_basin_mean_tech = groupsummary(mydata,["pottype","planttype","yrs","ssp","corner"],"mean","PnetAlls");

coe_subbasin_median_tech = groupsummary(mydata,["subname","pottype","planttype","yrs","ssp","corner"],"median","COEAlls");
coe_subbasin_mean_tech = groupsummary(mydata,["subname","pottype","planttype","yrs","ssp","corner"],"mean","COEAlls");

gwh_subbasin_median_tech = groupsummary(mydata,["subname","pottype","planttype","yrs","ssp","corner"],"median","PnetAlls");
gwh_subbasin_mean_tech = groupsummary(mydata,["subname","pottype","planttype","yrs","ssp","corner"],"mean","PnetAlls");
disp("Evaluated medians")

%% Plot median GWhs and COEs w color and symbol types distinguished - basin scale
title2=["GWhs" "COEs"];
figure
%idata=coe_basin_median;pp=2;
idata=gwh_basin_mean;pp=1;
gscatter(repmat(1:25,1, height(idata)/25,1),idata.median_COEAlls,{idata.pottype,idata.planttype},[cl_RP; cl_DP],"^^**")
set(gca,'XTick',1:25,'XTickLabel',["Historical" rcpcornernames rcpcornernames])
xline([1.5 13.5],"Color",mygraylines)
ylabel('Project median','FontWeight','bold')
grid on
grid minor
xlim([0 26])
applymyplotformat(title2(pp))
l=legend([" " " " planttypes{:}],'NumColumns',2);
l.Title.String=strcat(pottypes2{:});

%% Plot MEAN GWhs or COEs w color and symbol types distinguished - basin scale
title2=strcat("MEAN",["GWhs" "COEs"]);
figure
%idata=coe_basin_median;pp=2;
idata=gwh_basin_mean;pp=1;
gscatter(repmat(1:25,1, height(idata)/25,1),idata.mean_PnetAlls,{idata.pottype,idata.planttype},[cl_RP; cl_DP],"^^**")
set(gca,'XTick',1:25,'XTickLabel',["Historical" rcpcornernames rcpcornernames])
xline([1.5 13.5],"Color",mygraylines)
ylabel('Project MEAN','FontWeight','bold')
grid on
grid minor
xlim([0 26])
applymyplotformat(title2(pp))
l=legend([" " " " planttypes{:}],'NumColumns',2);
l.Title.String=strcat(pottypes2{:});

%% Plot SUBBASIN median GWhs  w color and symbol types distinguished -  All in one
title2=["GWhs" "COEs"];
selrows=gwh_subbasin_median.pottype=="Sustainable";
figure
%idata=coe_subbasin_median;pp=2;
idata=gwh_subbasin_median(selrows,:);pp=1;
%gscatter(repmat(1:25,1, height(idata)/25,1),idata.median_PnetAlls,{idata.subname, idata.pottype,idata.planttype},repelem(cmap8_wong_gray,4,1),"^so*")
gscatter(repmat(1:25,1, height(idata)/25,1),idata.median_PnetAlls,{idata.subname, idata.pottype,idata.planttype},[cl_RP; cl_DP],repelem('^so*.+xd',1,2))

set(gca,'XTick',1:25,'XTickLabel',["Historical" rcpcornernames rcpcornernames])
xline([1.5 13.5],"Color",mygraylines)
ylabel('Project median','FontWeight','bold')
grid on
grid minor
xlim([0 26])
applymyplotformat(title2(pp))
l=legend('NumColumns',2);
% l.Title.String=strcat(pottypes2{:});
set(gca,'YScale','log')

%% Plot SUBBASIN median GWhs  w color and symbol types distinguished -  All separate
title2=["GWhs" "COEs"];
ss=1;
figure
for sb=basindata.basinnames'
    selrows=gwh_subbasin_median.subname==sb; %gwh_subbasin_median.pottype=="Sustainable" &

    subplot(2,4,ss)%idata=coe_subbasin_median;pp=2;
    idata=gwh_subbasin_median(selrows,:);pp=1;
    %gscatter(repmat(1:25,1, height(idata)/25,1),idata.median_PnetAlls,{idata.subname, idata.pottype,idata.planttype},repelem(cmap8_wong_gray,4,1),"^so*")
    gscatter(repmat(1:25,1, height(idata)/25,1),idata.median_PnetAlls,{idata.pottype,idata.planttype},[cl_RP; cl_DP],'^ovs')
    set(gca,'XTick',1:25,'XTickLabel',["Historical" rcpcornernames rcpcornernames])
    xline([1.5 13.5],"Color",mygraylines)
    ylabel('Project median','FontWeight','bold')
    grid on
    grid minor
    xlim([0 26])
    applymyplotformat(sb)
    % l.Title.String=strcat(pottypes2{:});
    set(gca,'YScale','log')
    ss=ss+1;
    legend off
end
l=legend('NumColumns',3);
sgtitle(title2(pp))


%% Plot basin MEAN (magnitude and %) COEs and GWh in one plot w color and symbols
selfincols=[3,4];
selsustcol=[1,2];

figure;
colororder([cl_RP; cl_DP])
%GWH
basin_datain=gwh_basin_mean.mean_PnetAlls;
tmp_finsust=reshape(basin_datain,25,4);
tmp_finsust_prctchange=(tmp_finsust-tmp_finsust(1,:))./tmp_finsust(1,:)*100;
subplot(2,2,1)
hold all
scatter(1:25,tmp_finsust(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,tmp_finsust(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
set(gca,'XTick',1:25,'XTickLabel',"")
xline([1.5 13.5],"Color",mygraylines)
ylabel('Mean','FontWeight','bold')
xlim([0.5 25.5])
applymyplotformat("GWh - magnitude")
%grid minor
l=legend([" " " " planttypes{:}],'NumColumns',2,'Location','northoutside');
l.Title.String=strjoin(pottypes2(:),"  ");
set(gca,'YScale','log')%,'Ylim',[1, 200])
%GWH change
subplot(2,2,2)
hold all
scatter(1:25,tmp_finsust_prctchange(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,tmp_finsust_prctchange(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
set(gca,'XTick',1:25,'XTickLabel',"")
xline([1.5 13.5],"Color",mygraylines)
ylabel('% change from historical','FontWeight','bold')
yline(0,'Color','red')
xlim([1.5 25.5])
applymyplotformat("GWh - change")
%grid minor
xticklabels([])
text(6+[1.5,13.5],70*[1 1],tf_full,'HorizontalAlignment','center','FontAngle','italic')

%COE
basin_datain=coe_basin_mean.mean_COEAlls;
tmp_finsust=reshape(basin_datain,25,4);
tmp_finsust_prctchange=(tmp_finsust-tmp_finsust(1,:))./tmp_finsust(1,:)*100;
subplot(2,2,3)
hold all
scatter(1:25,tmp_finsust(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)%,'filled'
scatter(1:25,tmp_finsust(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)%,'filled'
set(gca,'XTick',1:25,'XTickLabel',["Historical" rcpcornernames rcpcornernames])%,'YColor', extraaxis)
xline([1.5 13.5],"Color",mygraylines)
ylabel('Mean','FontWeight','bold')
xlim([0.5 25.5])
applymyplotformat("COE - magnitude")
%grid minor
% COE change
subplot(2,2,4)
hold all
scatter(1:25,tmp_finsust_prctchange(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,tmp_finsust_prctchange(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
set(gca,'XTick',1:25,'XTickLabel',["Historical" rcpcornernames rcpcornernames])%,'YColor', extraaxis)
xline([1.5 13.5],"Color",mygraylines)
ylabel('% change from historical','FontWeight','bold')
yline(0,'Color','red')
xlim([1.5 25.5])
applymyplotformat("COE - change")
%grid minor

%% GOOD: Plot basin MEDIAN (magnitude and %) COEs and GWh in one plot w color and symbols
selfincols=[3,4];
selsustcol=[1,2];

figure;
colororder([cl_RP; cl_DP])
%GWH
basin_datain=gwh_basin_median.median_PnetAlls;
tmp_finsust=reshape(basin_datain,25,4);
tmp_finsust_prctchange=(tmp_finsust-tmp_finsust(1,:))./tmp_finsust(1,:)*100;
subplot(2,2,1)
hold all
scatter(1:25,tmp_finsust(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,tmp_finsust(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
set(gca,'XTick',1:25,'XTickLabel',"")
xline([1.5 13.5],"Color",mygraylines)
ylabel('Median','FontWeight','bold')
xlim([0.5 25.5])
applymyplotformat("GWh - magnitude")
%grid minor
l=legend([" " " " planttypes{:}],'NumColumns',2,'Location','northoutside');
if ff==1
    l.Title.String=strjoin(pottypes2(:),"  ");
    set(gca,'YScale','log','Ylim',[1, 200])

else
    l.Title.String=strjoin(pottypes3(2:3),"  ");
end
subplot(2,2,2)
hold all
scatter(1:25,tmp_finsust_prctchange(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,tmp_finsust_prctchange(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
set(gca,'XTick',1:25,'XTickLabel',"")
xline([1.5 13.5],"Color",mygraylines)
ylabel('% change from historical','FontWeight','bold')
yline(0,'Color','red')
xlim([1.5 25.5])
applymyplotformat("GWh - change")
%grid minor
xticklabels([])
topy=ylim;
text(6+[1.5,13.5],topy(2)*[1 1],tf_full,'HorizontalAlignment','center','FontAngle','italic')

%COE
basin_datain=coe_basin_median.median_COEAlls;
tmp_finsust=reshape(basin_datain,25,4);
tmp_finsust_prctchange=(tmp_finsust-tmp_finsust(1,:))./tmp_finsust(1,:)*100;
subplot(2,2,3)
hold all
scatter(1:25,tmp_finsust(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)%,'filled'
scatter(1:25,tmp_finsust(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)%,'filled'
set(gca,'XTick',1:25,'XTickLabel',["Historical" rcpcornernames rcpcornernames])%,'YColor', extraaxis)
xline([1.5 13.5],"Color",mygraylines)
ylabel('Median','FontWeight','bold')
xlim([0.5 25.5])
applymyplotformat("COE - magnitude")
%grid minor

subplot(2,2,4)
hold all
scatter(1:25,tmp_finsust_prctchange(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,tmp_finsust_prctchange(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
set(gca,'XTick',1:25,'XTickLabel',["Historical" rcpcornernames rcpcornernames])%,'YColor', extraaxis)
xline([1.5 13.5],"Color",mygraylines)
ylabel('% change from historical','FontWeight','bold')
yline(0,'Color','red')
xlim([1.5 25.5])
applymyplotformat("COE - change")
%grid minor


%% GOOD: Plot basin and subbasins MEDIAN (magnitude and %) GWh only w color and symbols - TECH, FIN AND SUST
ss=1; %subplot starting index
ff=3; % fig idx
%%Reshape data so easier to split RP and DP and calculate total and % change
%set basin data - tech
seltechcols=[1,2]; %RP and DP
basin_datain=gwh_basin_median_tech.median_PnetAlls;
tmp_tech=reshape(basin_datain,25,[]);
basin_tech_prctchange=(tmp_tech-tmp_tech(1,:))./tmp_tech(1,:)*100;

%set basin data - fin,sust
selfincols=[3,4];
selsustcol=[1,2];
basin_datain=gwh_basin_median.median_PnetAlls;
tmp_finsust=reshape(basin_datain,25,4);
basin_finsust_prctchange=(tmp_finsust-tmp_finsust(1,:))./tmp_finsust(1,:)*100;

%set subbasin data - tech, fin, sust
subbasin_datain_finsust=gwh_subbasin_median;
subbasin_datain_tech=gwh_subbasin_median_tech;

%GWH
figure(ff*100);clf;
colororder([cl_RP; cl_DP])
subplot(3,3,ss)
hold all
scatter(1:25,tmp_tech(:, seltechcols),30,'+','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,tmp_finsust(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,tmp_finsust(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
set(gca,'XTick',1:25,'XTickLabel',"",...
    'YScale','log')%,'Ylim',[1, 1000])
xline([1.5 13.5],"Color",mygraylines)
ylabel('Median','FontWeight','bold')
xlim([0.5 25.5])
applymyplotformat("GWh - magnitude ALL INDUS")
topy=ylim;
text(6+[1.5,13.5],topy(2)*[1 1],tf_full,'HorizontalAlignment','center','FontAngle','italic')
l=legend([" " " " " " " " planttypes{:}],'NumColumns',3,'Location','northoutside');
l.Title.String=strjoin(pottypes3,"  ");

%GWH change
figure(ff*200);clf;
colororder([cl_RP; cl_DP])
subplot(3,3,ss)
hold all
scatter(1:25,basin_tech_prctchange(:, seltechcols),30,'+','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,basin_finsust_prctchange(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,basin_finsust_prctchange(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
set(gca,'XTick',1:25,'XTickLabel',"")
xline([1.5 13.5],"Color",mygraylines)
ylabel('% change from historical','FontWeight','bold')
yline(0,'Color','red')
xlim([1.5 25.5]) % skip historical
applymyplotformat("GWh - change ALL INDUS")
%grid minor
xticklabels([])
topy=ylim;
text(6+[1.5,13.5],topy(2)*[1 1],tf_full,'HorizontalAlignment','center','FontAngle','italic')
l=legend([" " " " planttypes{:}],'NumColumns',2,'Location','northoutside');
if ff==1
    l.Title.String=strjoin(pottypes2(:),"  ");
else
    l=legend([" " " " " " " " planttypes{:}],'NumColumns',3,'Location','northoutside');
    l.Title.String=strjoin(pottypes3(1:3),"  ");
end

% Add subbasin MEDIAN % change GWh in on plot w color and symbols
ss=2;
for sb=basindata.basinnames'
    % for tech
    selrows=subbasin_datain_tech.subname==sb; %gwh_subbasin_median.pottype=="Sustainable" &
    %Calculate total and % change
    subbasin_datain2=subbasin_datain_tech.median_PnetAlls(selrows);
    tmp_tech=reshape(subbasin_datain2,25,[]);
    tmp_tech_prctchange=(tmp_tech-tmp_tech(1,:))./tmp_tech(1,:)*100;

    %for fin and sust
    selrows=subbasin_datain_finsust.subname==sb; %gwh_subbasin_median.pottype=="Sustainable" &

    %Calculate total and % change
    tmpsubbasin=subbasin_datain_finsust.median_PnetAlls(selrows);
    tmp_finsust=reshape(tmpsubbasin,25,[]);
    tmp_finsust_prctchange=(tmp_finsust-tmp_finsust(1,:))./tmp_finsust(1,:)*100;

    %GWH
    figure(ff*100)
    subplot(3,3,ss)
    hold all
    scatter(1:25,tmp_tech(:, seltechcols),30,'+','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    scatter(1:25,tmp_finsust(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    scatter(1:25,tmp_finsust(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    xline([1.5 13.5],"Color",mygraylines)
    ylabel('Median','FontWeight','bold')
    xlim([0.5 25.5])
    applymyplotformat(strcat("GWh - magnitude for ",sb))
    %grid minor
    set(gca,'YScale','log')
    if ff==1
        ylim([1, 1000])
    else
    end
    if ss>6
        set(gca,'XTick',1:25,'XTickLabel',["Historical" rcpcornernames rcpcornernames])%,'YColor', extraaxis)
    else
        set(gca,'XTick',1:25,'XTickLabel',"")
    end

    %GWH change
    figure(ff*200)
    subplot(3,3,ss)
    hold all
    scatter(1:25,tmp_tech_prctchange(:, seltechcols),30,'+','MarkerFaceAlpha',myalpha,"LineWidth",1.25)

    scatter(1:25,tmp_finsust_prctchange(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    scatter(1:25,tmp_finsust_prctchange(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    xline([1.5 13.5],"Color",mygraylines)
    ylabel('% change from historical','FontWeight','bold')
    yline(0,'Color','red')
    xlim([1.5 25.5])
    %ylim([-50 121])
    %ylim([-50 352])
    applymyplotformat(strcat("GWh - change for ", sb) )
    %grid minor
    xticklabels([])
    if ss>6
        set(gca,'XTick',1:25,'XTickLabel',["Historical" rcpcornernames rcpcornernames])%,'YColor', extraaxis)
    else
        set(gca,'XTick',1:25,'XTickLabel',"")
    end
    ss=ss+1;
end

%% GOOD: Plot basin and subbasins MEDIAN (magnitude and %) COE only w color and symbols - TECH, FIN AND SUST
ss=1; %subplot starting index
ff=3; % fig idx
%%Reshape data so easier to split RP and DP and calculate total and % change
%set basin data
seltechcols=[1,2]; %RP and DP
basin_datain=coe_basin_median_tech.median_COEAlls;
tmp_tech=reshape(basin_datain,25,[]);
basin_tech_prctchange=(tmp_tech-tmp_tech(1,:))./tmp_tech(1,:)*100;

selfincols=[3,4];
selsustcol=[1,2];
basin_datain=coe_basin_median.median_COEAlls;
tmp_finsust=reshape(basin_datain,25,4);
basin_finsust_prctchange=(tmp_finsust-tmp_finsust(1,:))./tmp_finsust(1,:)*100;

%set subbasin data
subbasin_datain_finsust=coe_subbasin_median;
subbasin_datain_tech=coe_subbasin_median_tech;

%COE
figure(ff*100);clf;
colororder([cl_RP; cl_DP])
subplot(3,3,ss)
hold all
scatter(1:25,tmp_tech(:, seltechcols),30,'+','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,tmp_finsust(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,tmp_finsust(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
set(gca,'XTick',1:25,'XTickLabel',"")
%,'YScale','log''Ylim',[1, 1000])
xline([1.5 13.5],"Color",mygraylines)
ylabel('Median','FontWeight','bold')
xlim([0.5 25.5])
applymyplotformat("COE - magnitude ALL INDUS")
topy=ylim;
text(6+[1.5,13.5],topy(2)*[1 1],tf_full,'HorizontalAlignment','center','FontAngle','italic')
l=legend([" " " " " " " " planttypes{:}],'NumColumns',3,'Location','northoutside');
l.Title.String=strjoin(pottypes3,"  ");

%COE change
figure(ff*200);clf;
colororder([cl_RP; cl_DP])
subplot(3,3,ss)
hold all
scatter(1:25,basin_tech_prctchange(:, seltechcols),30,'+','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,basin_finsust_prctchange(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,basin_finsust_prctchange(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
set(gca,'XTick',1:25,'XTickLabel',"")
xline([1.5 13.5],"Color",mygraylines)
ylabel('% change from historical','FontWeight','bold')
yline(0,'Color','red')
xlim([1.5 25.5]) % skip historical
applymyplotformat("GWh - change ALL INDUS")
%grid minor
xticklabels([])
topy=ylim;
text(6+[1.5,13.5],topy(2)*[1 1],tf_full,'HorizontalAlignment','center','FontAngle','italic')
l=legend([" " " " planttypes{:}],'NumColumns',2,'Location','northoutside');
if ff==1
    l.Title.String=strjoin(pottypes2(:),"  ");
else
    l=legend([" " " " " " " " planttypes{:}],'NumColumns',3,'Location','northoutside');
    l.Title.String=strjoin(pottypes3(1:3),"  ");
end

% Add subbasin MEDIAN % change GWh in on plot w color and symbols
ss=2;
for sb=basindata.basinnames'
    % for tech
    selrows=subbasin_datain_tech.subname==sb; %gwh_subbasin_median.pottype=="Sustainable" &
    %Calculate total and % change
    subbasin_datain2=subbasin_datain_tech.median_COEAlls(selrows);
    tmp_tech=reshape(subbasin_datain2,25,[]);
    tmp_tech_prctchange=(tmp_tech-tmp_tech(1,:))./tmp_tech(1,:)*100;

    %for fin and sust
    selrows=subbasin_datain_finsust.subname==sb; %gwh_subbasin_median.pottype=="Sustainable" &

    %Calculate total and % change
    tmpsubbasin=subbasin_datain_finsust.median_COEAlls(selrows);
    tmp_finsust=reshape(tmpsubbasin,25,[]);
    tmp_finsust_prctchange=(tmp_finsust-tmp_finsust(1,:))./tmp_finsust(1,:)*100;

    %GWH
    figure(ff*100)
    subplot(3,3,ss)
    hold all
    scatter(1:25,tmp_tech(:, seltechcols),30,'+','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    scatter(1:25,tmp_finsust(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    scatter(1:25,tmp_finsust(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    xline([1.5 13.5],"Color",mygraylines)
    ylabel('Median','FontWeight','bold')
    xlim([0.5 25.5])
    applymyplotformat(strcat("COE - magnitude for ",sb))
    %grid minor
    %set(gca,'YScale','log')
    if ff==1
        ylim([1, 1000])
    else
    end
    if ss>6
        set(gca,'XTick',1:25,'XTickLabel',["Historical" rcpcornernames rcpcornernames])%,'YColor', extraaxis)
    else
        set(gca,'XTick',1:25,'XTickLabel',"")
    end

    %GWH change
    figure(ff*200)
    subplot(3,3,ss)
    hold all
    scatter(1:25,tmp_tech_prctchange(:, seltechcols),30,'+','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    scatter(1:25,tmp_finsust_prctchange(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    scatter(1:25,tmp_finsust_prctchange(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    xline([1.5 13.5],"Color",mygraylines)
    ylabel('% change from historical','FontWeight','bold')
    yline(0,'Color','red')
    xlim([1.5 25.5])
    %ylim([-50 121])
    %ylim([-50 352])
    applymyplotformat(strcat("COE - change for ", sb) )
    %grid minor
    xticklabels([])
    if ss>6
        set(gca,'XTick',1:25,'XTickLabel',["Historical" rcpcornernames rcpcornernames])%,'YColor', extraaxis)
    else
        set(gca,'XTick',1:25,'XTickLabel',"")
    end
    ss=ss+1;
end



%% Write min max for each RCP % subbasin scale to sheet
oxlsfile=fullfile(rootoffut,"RangeinDPRP_COE_GWh.xlsx");
basingrp=["pottype","planttype","yrs","ssp"];
basin_all=[coe_basin_median; coe_basin_median_tech];
basin_all.median_PnetAlls=[gwh_basin_median.median_PnetAlls; gwh_basin_median_tech.median_PnetAlls];
summary_basin = groupsummary(basin_all,basingrp,["min","max"],["median_COEAlls" "median_PnetAlls"]);
summary_basin_norcp = flipud(groupsummary(basin_all,basingrp([2,1,3]),["min","max"],["median_COEAlls" "median_PnetAlls"]));

subbasin_all=[coe_subbasin_median; coe_subbasin_median_tech];
subbasin_all.median_PnetAlls=[gwh_subbasin_median.median_PnetAlls; gwh_subbasin_median_tech.median_PnetAlls];
subbasingrp=["pottype","subname","planttype","yrs","ssp"];
summary_basin = groupsummary(subbasin_all,subbasingrp,["min","max"],["median_COEAlls" "median_PnetAlls"]);

writetable(summary_basin,oxlsfile,'Sheet','Basin_vals')
writetable(summary_basin_norcp,oxlsfile,'Sheet','Basin_vals_norcp')

writetable(summary_subbasin,oxlsfile,'Sheet','SubBasin_vals')

disp("Written min-max of medians to excel")


