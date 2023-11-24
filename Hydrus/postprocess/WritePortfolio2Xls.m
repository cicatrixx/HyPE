%% WritePortfolio2Xls
% Write main outputs relevant for CAS to excel in previous template
% Also write historical portfolio for Arthur.

clc
clearvars
close all
addpath('Hydrus\postprocess') % Add folder path
run('myVarNames.m')
printsummary=1;
CASxlspath=fullfile(rootf,'output','ForCAS','CAS_01Mar2022'); if ~isfolder(CASxlspath); mkdir(CASxlspath); end
xlsVarnames={'PID', 'Latss', 'Lonss', 'SysIDs', 'COEAlls_[USD2010/kWh]', 'PnetAlls_[GWh/yr]', 'BasinID'};

%% Write the tech2sust portfolios
load('G:\SurfDrive\HPmodel\output\HistRuns_Figs_Analysis\Figs_trial\MainScenarios.mat', 'pcsoutT2S','runnamesT2S')
runname_prefix='R103_';
selpot="Tech2Sust";
for fid=2:3 % only load the sust water and sust land portfolios
    seldata=pcsoutT2S{fid};
    selecon= seldata.COEAlls<=costlim & ~isnan(seldata.COEAlls);
    outdata=table(seldata.PID(selecon),seldata.latss(selecon),seldata.lonss(selecon),...
        seldata.SysIDs(selecon),seldata.COEAlls(selecon),seldata.PnetAlls(selecon),...
        seldata.subbasinIDs(selecon),'VariableNames',xlsVarnames);  % exported to excel

    % display sum for quick check
    if printsummary
        disp(runnamesT2S{fid})
        fprintf('%s Potential: %0.2f TWh\n',selpot, nansum(seldata.PnetAlls(seldata.COEAlls<=0.1))/1000)
        fprintf('# of selected projects: %0.0f out of %0.0f\n',sum(~isnan(seldata.PnetAlls(seldata.COEAlls<=0.1))),length(seldata.PnetAlls))
        fprintf('# of selected River type projects: %0.0f with %0.2f GWh\n',sum(seldata.SysIDs==2 & seldata.COEAlls<=0.1),sum(seldata.PnetAlls(seldata.SysIDs==2 & seldata.COEAlls<=0.1)))
        fprintf('# of selected Diversion type projects: %0.0f with %0.2f GWh\n',sum(seldata.SysIDs==1 & seldata.COEAlls<=0.1),sum(seldata.PnetAlls(seldata.SysIDs==1 & seldata.COEAlls<=0.1)))
    end
    % write to excel
    oxlsfile=fullfile(CASxlspath,['R103_Tech2Sust_Full_Mixed.xlsx']);
    sheetname=extractBetween(runnamesT2S{fid},runname_prefix,"_Full");
    writetable(outdata,oxlsfile,'Sheet',sheetname{:});
end
disp("Written historical Tech2Sust scenario portfolio totals to excel")

%% Write the future scenario totals
run('myVarNames_Fut.m')
ofname='Compile_FutScenario_Totals00';
matfile_summary=fullfile(rootf,[ofname '.mat']);
xlsfile_summary=fullfile(CASxlspath,[ofname '.xlsx']);

% Load and clean up previous compilation for CAS use. Remove unnecessary
% rows and columns, use more intuitive names.
load(fullfile(rootof,[ofname '.mat']), 'tot_allscen')
tot_allscen.rcp=renamecats(tot_allscen.ssp,sspshort,rcpnames);
tot_allscen.futterms=tot_allscen.yrs;

reordervar= ["Runname","pottype","futterms", "rcp","corner","gcm","TWh/yr of Diversion canal","TWh/yr of River power","Total TWh/yr",...
    "Number of Diversion canal", "Number of River power","Total number", 'Prct nplants costing <=0.10'...
    'Historical TWh/yr', 'Historical Num Proj', 'Historical Prct nplants costing <=0.10',    ];
tot_allscen=tot_allscen(:,reordervar);
save(matfile_summary,'tot_allscen')

notSustTech=tot_allscen.pottype ~= 'Sust-Tech';
writetable(tot_allscen(notSustTech,1:12),xlsfile_summary,'Sheet',"TotalsCompiled");
disp("Written all future scenario portfolio totals to excel")

%% Write the future portfolios - each future term in separate folder
load(fullfile(rootof,'MainFutScenarios.mat'), 'pcsout', 'runnames')
selterm=1;
namecnt=1;
rcpcornerfnames=strrep(strrep(strrep(rcpcornernames,".","_"),": ","_")," ","");
% First 2 scenarios are historical so can skip
% Loop through tech/sust for each scenario
for fid=3:2:length(pcsout)
    % Create separate folder for each term
    if fid==27; namecnt=1; selterm=2; end
    if fid==3 ||fid==27
        mkdir(fullfile(CASxlspath,termname(selterm)))
    end

    %currentrun=newfnames{fidx};
    currentrun=rcpcornerfnames{namecnt};
    namecnt=namecnt+1;  % for each scenario

    % Eval tech and econ data
    seldata=pcsout{fid};
    seltech=~isnan(seldata.COEAlls);
    selecon= seldata.COEAlls<=costlim & seltech;

    tech_outdata_ss=table(seldata.PID(seltech),seldata.latss(seltech),seldata.lonss(seltech),...
        seldata.SysIDs(seltech),seldata.COEAlls(seltech),seldata.PnetAlls(seltech),...
        seldata.subbasinIDs(seltech),'VariableNames',xlsVarnames);

    fin_outdata_ss=table(seldata.PID(selecon),seldata.latss(selecon),seldata.lonss(selecon),...
        seldata.SysIDs(selecon),seldata.COEAlls(selecon),seldata.PnetAlls(selecon),...
        seldata.subbasinIDs(selecon),'VariableNames',xlsVarnames);  % exported to excel

    % Eval sust data
    selsustdata=pcsout{fid+1};
    selsust= selsustdata.COEAlls<=costlim & ~isnan(selsustdata.COEAlls);

    sust_outdata_ss=table(selsustdata.PID(selsust),selsustdata.latss(selsust),selsustdata.lonss(selsust),...
        selsustdata.SysIDs(selsust),selsustdata.COEAlls(selsust),selsustdata.PnetAlls(selsust),...
        selsustdata.subbasinIDs(selsust),'VariableNames',xlsVarnames);  % exported to excel
    % write to excel
    oxlsfile=fullfile(CASxlspath,termname(selterm),[currentrun '.xlsx']);
    writetable(tech_outdata_ss,oxlsfile,'Sheet','Tech_RiskAverse')
    writetable(fin_outdata_ss,oxlsfile,'Sheet','Fin_RiskAverse')
    writetable(sust_outdata_ss,oxlsfile,'Sheet','Sust_RiskAverse')

end
disp("Written all future scenario portfolios to excel")

%% Write historical portfolio - only mixed tech,fin,sust



%% Write the protected area layers
load('G:\SurfDrive\HPmodel\data\ASIA\Basin_UIB\PantpeBasin_103.mat', 'Heritages_cultural', 'Heritages_natural')
%buffers
cellsz_m=500; 
buffer_cult_km = 1;         %Buffer for cultural heritages in km
buffer_nat_km = 2;          %Buffer for natural heritages in km
ProtArea = createBuffer(Heritages_natural, buffer_nat_km*1e3/cellsz_m)+ ...
    createBuffer(Heritages_cultural, buffer_cult_km*1e3/cellsz_m);

savemat2Pantpetiff(fullfile(rootf, 'data', 'ASIA', 'Basin_UIB',"protectedareas_wbuffer.tif"), ProtArea)

%% Manually ran parts of runhydrus to recreate and save geohazard layer for CAS visuals
% savemat2Pantpetiff(fullfile(rootf, 'data', 'ASIA', 'Basin_UIB',"geohazard_riskaverse_skippedhighriskarea.tif"), geohazard_skip)
% savemat2Pantpetiff(fullfile(rootf, 'data', 'ASIA', 'Basin_UIB',"geohazard_riskaverse_addedcost.tif"), geohazard_hazard_rates)


