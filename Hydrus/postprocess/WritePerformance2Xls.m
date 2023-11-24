% WritePerformance2XLs.m
% Created By    : Sanita Dhaubanjar
% Created For	: SustainIndus project WP2 Paper publication
%=========================
% Code for saving annual hydropower generation time-series from mat file to
% excel for Arthur. Only tech and sust as did not save fin

clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

%subplottight = @(m,n,p) subtightplot (m, n, p, [0.08 0.025],[.08 .08],[.08 .01]); %[vgap hgap], marg_h -[low up],marg_w -[l r]

%% Load HP performance TS
matname=fullfile(rootof,'FutScenarios_HPgen_basintot03'); %w/o extension as it is added later
selnewpots=[1 4];
pot2=pottypes3_short([1,3]);
for pot=1:2 % 1=tech, 2=sust
    HP_TS{pot}=load(strcat(matname,"_",newpots4{selnewpots(pot)},".mat"),'hist_annualTS_basin','fut_annualTS_basin',...
        "hist_monTS_basin","fut_annualTS_basin", "fut_monTS_basin", "cl_Qnames", "cl_pfnames" ); % most of these have data in the form of (pf, timestep, qf)
    % Convert to twh
    HP_TS{pot}.hist_annualTS_basin=HP_TS{pot}.hist_annualTS_basin/1000;
    HP_TS{pot}.fut_annualTS_basin=HP_TS{pot}.fut_annualTS_basin/1000;
    HP_TS{pot}.cl_pfnames(1)="Historical";

    fprintf("Loaded HP annual TS for %s pot \n", newpots4{selnewpots(pot)})
end

futmon_ts1=futmon_ts1';
futmon_ts2=futmon_ts2';

%% FINAL: TS of each pf under histQ and its own futQ
for pot=1%:2
    figure
    for pfidx=1:25
        nexttile
        %pf under histQ
        l1(1)= plot(histyrs_ts,HP_TS{pot}.hist_annualTS_basin(pfidx,:,1)','k-','DisplayName',"Histpf",'LineWidth',1.5);
        hold all
        % mid future
        %for qfidx=2:13

        if pfidx~=1
            %pf under its own futQ except for hist pf
            qfidx=pfidx;
            l1(qfidx)= plot(fyrs_ts1,HP_TS{pot}.fut_annualTS_basin(pfidx,:,qfidx)','DisplayName',"Histpf",'LineWidth',1.5);
        end
        applymyplotformat(strcat("PF: ", HP_TS{pot}.cl_pfnames(pfidx)))
        ylabel("Hydropower (TWh/yr)")
        % midfar future labels
        xline([2036, 2066],'LineStyle',':')
        topy=ylim;
        text([2050 2080],[1 1]*floor(topy(2)-50),tname,'HorizontalAlignment','center','FontWeight','bold')
    end
end

lgd=legend(["Historical", "Future"],'Orientation','horizontal');
lgd.Title.String="Porfolios";

%% FINAL: Write to excel TS of each pf under histQ and its own futQ
ALxlspath=fullfile(rootf,'output','ForAL_11Sep2023');
if ~isfolder(ALxlspath); mkdir(ALxlspath); end

for pot=1:2
    compileHistTS=HP_TS{pot}.hist_annualTS_basin(:,:,1);
    compileFutTS=zeros(25,30);

    for pfidx=1:25
        if pfidx~=1
            %pf under its own futQ except for hist pf
            qfidx=pfidx;
            compileFutTS(pfidx,:)=HP_TS{pot}.fut_annualTS_basin(pfidx,:,qfidx);
        end
    end
    %% Write to excel
    outdata=table(HP_TS{pot}.cl_pfnames,compileFutTS,'VariableNames',{'Portfolio','Future TS 2036-2065 or 2066-2095 '});
    oxlsfile=fullfile(ALxlspath,'FutmonTS_HPgen_basintot03.xlsx');
    writetable(outdata,oxlsfile,'Sheet',strcat(pot2{pot},'_TWh_yr'));
end

disp("Written future monthly TS to excel")