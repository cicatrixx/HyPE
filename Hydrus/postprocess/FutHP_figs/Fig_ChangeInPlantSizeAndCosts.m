% Fig_ChangeInPlantSizeAndCosts.m
% Created By    : Sanita Dhaubanjar
% Created For	: SustainIndus project WP2 Paper publication
%=========================
% Code for loading pre-processed output files from Hydrus model to create
% final plots used for publication
% Plot change in median plant size and plant cost for tech, fin, sust
% potential

clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

ofmatname=fullfile(rootof,"FutFigs_cl0","FutFigs_cl1","FigData",'MainFutScenarios.mat');
load(ofmatname, 'pcsout',...
    'runnames','basindata','catchments_cl')

%% Compile scenarios details
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
all_portfolios_compiled=table();
for scn=1:50
    scen_data_tmp=struct2table(rmfield(pcsout{scn},{'PID','co_arc','ro_arc','runname'}));
    %table(pcsout{scn}.PnetAlls,pcsout{scn}.COEAlls, pcsout{scn}.SysIDs,pcsout{scn}.subbasinIDs,...
    % 'VariableNames',{'GWh/yr','USD/kWh','PlantType','Subbasin'});

    %remove rows with missing data
    scen_data_tmp=rmmissing(scen_data_tmp);
    scen_data_tmp(:,'runname')={pcsout{scn}.runname};
    scen_data_tmp=[scen_data_tmp repmat(scen_details(scn,:),height(scen_data_tmp),1)];
    all_portfolios_compiled=[all_portfolios_compiled; scen_data_tmp];
end

all_portfolios_compiled.subname= categorical(all_portfolios_compiled.subbasinIDs,basindata.basinIDs,basindata.basinnames);
all_portfolios_compiled.planttype= categorical(all_portfolios_compiled.SysIDs,[2 1], planttypes);

disp("Created long table w HP data")
selfinpot=all_portfolios_compiled.COEAlls<=costlim;

%% Get medians and means for all or only financial project
onlyfinprjs=01;
if ~onlyfinprjs
    % this one would be technical pot
    coe_basin_median = groupsummary(all_portfolios_compiled,["pottype","planttype","yrs","ssp","corner"],"median","COEAlls");
    coe_basin_mean = groupsummary(all_portfolios_compiled,["pottype","planttype","yrs","ssp","corner"],"mean","COEAlls");

    gwh_basin_median = groupsummary(all_portfolios_compiled,["pottype","planttype","yrs","ssp","corner"],"median","PnetAlls");
    gwh_basin_mean = groupsummary(all_portfolios_compiled,["pottype","planttype","yrs","ssp","corner"],"mean","PnetAlls");
    %gwh_basin_mean = groupsummary(scen_data.PnetAlls,scen_data{:,["pottype","planttype","yrs","ssp","corner"]},"mean");


    coe_subbasin_median = groupsummary(all_portfolios_compiled,["subname","pottype","planttype","yrs","ssp","corner"],"median","COEAlls");
    coe_subbasin_mean = groupsummary(all_portfolios_compiled,["subname","pottype","planttype","yrs","ssp","corner"],"mean","COEAlls");

    gwh_subbasin_median = groupsummary(all_portfolios_compiled,["subname","pottype","planttype","yrs","ssp","corner"],"median","PnetAlls");
    gwh_subbasin_mean = groupsummary(all_portfolios_compiled,["subname","pottype","planttype","yrs","ssp","corner"],"mean","PnetAlls");
    ff=1; % idx for separating figs
else
    mydata=all_portfolios_compiled(selfinpot,:);
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
selrows=all_portfolios_compiled.pottype=="Technical";
mydata=all_portfolios_compiled(selrows,:);
coe_basin_median_tech = groupsummary(mydata,["pottype","planttype","yrs","ssp","corner"],"median","COEAlls");
coe_basin_mean_tech = groupsummary(mydata,["pottype","planttype","yrs","ssp","corner"],"mean","COEAlls");

gwh_basin_median_tech = groupsummary(mydata,["pottype","planttype","yrs","ssp","corner"],"median","PnetAlls");
gwh_basin_mean_tech = groupsummary(mydata,["pottype","planttype","yrs","ssp","corner"],"mean","PnetAlls");

coe_subbasin_median_tech = groupsummary(mydata,["subname","pottype","planttype","yrs","ssp","corner"],"median","COEAlls");
coe_subbasin_mean_tech = groupsummary(mydata,["subname","pottype","planttype","yrs","ssp","corner"],"mean","COEAlls");

gwh_subbasin_median_tech = groupsummary(mydata,["subname","pottype","planttype","yrs","ssp","corner"],"median","PnetAlls");
gwh_subbasin_mean_tech = groupsummary(mydata,["subname","pottype","planttype","yrs","ssp","corner"],"mean","PnetAlls");
disp("Evaluated medians")

%% GOOD: Plot basin and subbasins MEDIAN (magnitude and %) GWh only w color and symbols - TECH, FIN AND SUST
subplottight = @(m,n,p) subtightplot (m, n, p, [0.03 0.04],[.1 .09],[.08 .01]); %[vgap hgap], marg_h -[low up],marg_w -[l r]

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
paramname="% change in median plant size";
spidx=1; %subplot starting index

%GWH change in percentage
figure
colororder([cl_RP; cl_DP])
subplottight(3,3,spidx)
hold all
scatter(1:25,basin_tech_prctchange(:, seltechcols),30,'+','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,basin_finsust_prctchange(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,basin_finsust_prctchange(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
set(gca,'XTick',1:25,'XTickLabel',"")
xline([1.5 13.5],"Color",mygraylines)
%ylabel(paramname)
yline(0,'Color','red')
xlim([1.5 25.5]) % skip historical
applymyplotformat(strcat(charlbl(spidx),label_subbasin_basin(9)))
%grid minor
xticklabels([])
topy=ylim;
%text(6+[1.5,13.5],topy(2)*[1 1],tf_full,'HorizontalAlignment','center','FontAngle','italic')
text(6+[1.5,13.5],topy(2)*[1 1],tname,'HorizontalAlignment','center') %,'FontAngle','italic')

l=legend([" " " " " " " " planttypes{:}],'NumColumns',3,'Location','northoutside');
l.Title.String=strjoin(pottypes3(1:3),"  ");

% Add subbasin MEDIAN % change GWh in on plot w color and symbols
spidx=2;
for sb=basindata.basinnames'
    % Calculate total and % change for tech
    selrows=subbasin_datain_tech.subname==sb; %gwh_subbasin_median.pottype=="Sustainable" &
    subbasin_datain2=subbasin_datain_tech.median_PnetAlls(selrows);
    tmp_tech=reshape(subbasin_datain2,25,[]);
    tmp_tech_prctchange=(tmp_tech-tmp_tech(1,:))./tmp_tech(1,:)*100;

    %Calculate total and % change for fin and sust
    selrows=subbasin_datain_finsust.subname==sb; %gwh_subbasin_median.pottype=="Sustainable" &
    tmpsubbasin=subbasin_datain_finsust.median_PnetAlls(selrows);
    tmp_finsust=reshape(tmpsubbasin,25,[]);
    tmp_finsust_prctchange=(tmp_finsust-tmp_finsust(1,:))./tmp_finsust(1,:)*100;

    %GWH change
    subplottight(3,3,spidx)
    hold all
    scatter(1:25,tmp_tech_prctchange(:, seltechcols),30,'+','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    scatter(1:25,tmp_finsust_prctchange(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    scatter(1:25,tmp_finsust_prctchange(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    xline([1.5 13.5],"Color",mygraylines)

    yline(0,'Color','red')
    xlim([1.5 25.5])
    %ylim([-50 121])
    %ylim([-50 352])
    applymyplotformat(strcat(charlbl(spidx),sb))

    % ylabel in left middle plot
    if spidx==4
        ylabel(paramname)
    end

    %grid minor
    xticklabels([])
    if spidx>6        
        set(gca,'XTick',1:25,'XTickLabel',["Historical" cornernamesrep cornernamesrep] )%,'YColor', extraaxis)
        % add rcpnames below corners
        topy=ylim;
        text([3.5:4:25.5],topy(1)*ones(1,6),strrep([rcpnames rcpnames]," ",""),'HorizontalAlignment','center','fontweight','bold','FontSize',11);
    else
        set(gca,'XTick',1:25,'XTickLabel',"")
    end
    spidx=spidx+1;
end

%% GOOD: Plot basin and subbasins MEDIAN (magnitude and %) COE only w color and symbols - TECH, FIN AND SUST
subplottight = @(m,n,p) subtightplot (m, n, p, [0.03 0.04],[.1 .09],[.08 .01]); %[vgap hgap], marg_h -[low up],marg_w -[l r]

%%Reshape data so easier to split RP and DP and calculate total and % change
%set basin data - tech
seltechcols=[1,2]; %RP and DP
basin_datain=coe_basin_median_tech.median_COEAlls;
tmp_tech=reshape(basin_datain,25,[]);
basin_tech_prctchange=(tmp_tech-tmp_tech(1,:))./tmp_tech(1,:)*100;

%set basin data - fin,sust
selfincols=[3,4];
selsustcol=[1,2];
basin_datain=coe_basin_median.median_COEAlls;
tmp_finsust=reshape(basin_datain,25,4);
basin_finsust_prctchange=(tmp_finsust-tmp_finsust(1,:))./tmp_finsust(1,:)*100;

%set subbasin data - tech, fin, sust
subbasin_datain_finsust=coe_subbasin_median;
subbasin_datain_tech=coe_subbasin_median_tech;
paramname="% change in median plant cost";
spidx=1; %subplot starting index
ff=3; % fig idx

%coe change in percentage
figure
colororder([cl_RP; cl_DP])
subplottight(3,3,spidx)
hold all
scatter(1:25,basin_tech_prctchange(:, seltechcols),30,'+','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,basin_finsust_prctchange(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
scatter(1:25,basin_finsust_prctchange(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
set(gca,'XTick',1:25,'XTickLabel',"")
xline([1.5 13.5],"Color",mygraylines)
%ylabel(paramname)
yline(0,'Color','red')
xlim([1.5 25.5]) % skip historical
applymyplotformat(strcat(charlbl(spidx),label_subbasin_basin(9)))
%grid minor
xticklabels([])
topy=ylim;
%text(6+[1.5,13.5],topy(2)*[1 1],tf_full,'HorizontalAlignment','center','FontAngle','italic')
text(6+[1.5,13.5],topy(2)*[1 1],tname,'HorizontalAlignment','center') %,'FontAngle','italic')

l=legend([" " " " " " " " planttypes{:}],'NumColumns',3,'Location','northoutside');
l.Title.String=strjoin(pottypes3(1:3),"  ");

% Add subbasin MEDIAN % change GWh in on plot w color and symbols
spidx=2;
for sb=basindata.basinnames'
    % Calculate total and % change for tech
    selrows=subbasin_datain_tech.subname==sb; %gwh_subbasin_median.pottype=="Sustainable" &
    subbasin_datain2=subbasin_datain_tech.median_COEAlls(selrows);
    tmp_tech=reshape(subbasin_datain2,25,[]);
    tmp_tech_prctchange=(tmp_tech-tmp_tech(1,:))./tmp_tech(1,:)*100;

    %Calculate total and % change for fin and sust
    selrows=subbasin_datain_finsust.subname==sb; %gwh_subbasin_median.pottype=="Sustainable" &
    tmpsubbasin=subbasin_datain_finsust.median_COEAlls(selrows);
    tmp_finsust=reshape(tmpsubbasin,25,[]);
    tmp_finsust_prctchange=(tmp_finsust-tmp_finsust(1,:))./tmp_finsust(1,:)*100;

    %GWH change
    subplottight(3,3,spidx)
    hold all
    scatter(1:25,tmp_tech_prctchange(:, seltechcols),30,'+','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    scatter(1:25,tmp_finsust_prctchange(:, selfincols),30,'^','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    scatter(1:25,tmp_finsust_prctchange(:, selsustcol),30,'o','MarkerFaceAlpha',myalpha,"LineWidth",1.25)
    xline([1.5 13.5],"Color",mygraylines)

    yline(0,'Color','red')
    xlim([1.5 25.5])
    %ylim([-50 121])
    %ylim([-50 352])
    applymyplotformat(strcat(charlbl(spidx),sb))

    % ylabel in left middle plot
    if spidx==4
        ylabel(paramname)
    end

    %grid minor
    xticklabels([])
    if spidx>6
        set(gca,'XTick',1:25,'XTickLabel',["Historical" cornernamesrep cornernamesrep])%,'YColor', extraaxis)
        % add rcpnames below corners
        topy=ylim;
        text([3.5:4:25.5],topy(1)*ones(1,6),strrep([rcpnames rcpnames]," ",""),'HorizontalAlignment','center','fontweight','bold','FontSize',11);
    else
        set(gca,'XTick',1:25,'XTickLabel',"")
    end
    spidx=spidx+1;
end

%%

disp("***************************************EOF***************************************")

