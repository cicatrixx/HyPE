% Plot figures for delta change in total energy and number of projects with
% change in design discharge values. Figs show delta in terms of magnitude,
% or % of total potential or % of RP/DP potentials.
% saves delta change in magnitude and as % of total to excel

clc
clearvars
close all
addpath('Hydrus\postprocess') % Add folder path

%addpath(pwd,genpath('..\GitConnect\Hydrus\devFiles\'), ...
%    genpath('Hydrus\'))
run('myVarNames.m')

% Load data
load(fullfile(pwd,'\output\Compile_Scenario_Totals03.mat'),'tot_allscen',...
    'mainidx_end','SA_Sust_idx_end','selmixsust')
%  save(matfile_summary,'tot_allscen', 'SA_Others_delta','mainidx_end','SA_Sust_idx_end',...
%       'SA_Tech_idx_end', 'SA_Fin_idx_end','selmixsust')

% FINAL SA_SustQ: Plot circle matrix of total val circles SA_SustQ
RP_Qs =      {'Q30', 'Q40', 'Q50'};
largeDP_Qs = {'Q25', 'Q30', 'Q40'};
smallDP_Qs = {'Q70', 'Q80', 'Q90'};
simpleplots=01;
saveout2xl=10;
oxlfname='Compile_Scenario_Totals03';

%% SA_SustQ: Get delta change
sel_SA_Q_idx=mainidx_end+1:SA_Sust_idx_end(1);
SA_Q_varnames=extractBefore(tot_allscen{sel_SA_Q_idx,1},'-Full');

% delta = scen - default mixed sust scen
SA_Q_delta=tot_allscen{sel_SA_Q_idx, 2:3}-tot_allscen{selmixsust,2:3};
SA_Q_delta_num=tot_allscen{sel_SA_Q_idx, 5:6}-tot_allscen{selmixsust,5:6};

%% SA_SustQ: Plot simple bar for all delta in magnitude and % for SA Q scenarios to get suitable plot limits
if simpleplots
    figure;
    subplot(2,1,1)
    b1=bar(categorical(SA_Q_varnames,SA_Q_varnames),SA_Q_delta,'stacked','FaceAlpha',baralpha);
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end
    ylim([-21 13])
    ylabel('\Delta TWh/yr')


    yyaxis right
    scatter(categorical(SA_Q_varnames,SA_Q_varnames),SA_Q_delta_num(:,1),50,cl_RP*.6,'d','filled') %,'MarkerFaceAlpha',baralpha)
    hold on
    scatter(categorical(SA_Q_varnames,SA_Q_varnames),SA_Q_delta_num(:,1)+SA_Q_delta_num(:,2),50,cl_DP*.6,'d','filled') % ,'MarkerFaceAlpha',baralpha) %manually stack scatter
    ylim([-210 130])

    legend(tot_allscen.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end
    ax = gca;
    ax.YAxis(2).Color = 'k';
    grid on
    title("Delta change for SA Qs Sust")

    % SA_SustQ: Plot simple bar for delta % SA Q scenarios to get suitable plot limits
    subplot(2,1,2)
    b1=bar(categorical(SA_Q_varnames,SA_Q_varnames),SA_Q_delta_prct,'stacked','FaceAlpha',baralpha);
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end
    %ylim([-45 35])

    ylabel('\Delta in %')


    %yyaxis right
    hold on
    scatter(categorical(SA_Q_varnames,SA_Q_varnames),SA_Q_delta_num_prct(:,1),50,cl_RP*.6,'d','filled') %,'MarkerFaceAlpha',baralpha)
    scatter(categorical(SA_Q_varnames,SA_Q_varnames),SA_Q_delta_num_prct(:,1)+SA_Q_delta_num_prct(:,2),50,cl_DP*.6,'d','filled') % ,'MarkerFaceAlpha',baralpha) %manually stack scatter
    %ylim([-45 35])

    %legend(tot_allscen.Properties.VariableNames{[2:3,5:6]})%,'Location','northoutside','Orientation','horizontal')%)
    for k = 1:2
        b1(k).FaceColor = mycolors{k};
        b1(k).EdgeColor = mycolors{k};
    end
    %ax = gca;
    %ax.YAxis(2).Color = 'k';
    grid on
    title("Delta change in percentage for SA Qs Sust")
end

%% FINAL SA_SustQ: Plot bar matrix for delta SA Q scenarios
tt=1;
runcnt=1;
f=figure; %subplot(3,1,2)
for ii=1:numel(RP_Qs)
    hold all
    for jj=1:numel(largeDP_Qs)
        subtightplot (3, 3, tt, [0.025 0.08],[.05 .1],[.1 .06]); %[vgap hgap], marg_h,marg_w subplot(3,3,tt);
        selidx=runcnt:runcnt+2;
        b1=bar(categorical(SA_Q_varnames(selidx),SA_Q_varnames(selidx)),SA_Q_delta(selidx,:),'stacked','FaceAlpha',baralpha,'EdgeAlpha',0);

        for k = 1:2
            b1(k).FaceColor = mycolors{k};
            b1(k).EdgeColor = mycolors{k};
        end
        ylim([-21 14])
        if jj==1 % add left ylabel to left col
            ylabel('\Delta TWh/yr')
            text(-0.8,-5.5,sprintf('RP = %s',RP_Qs{ii}),'FontWeight','bold','Rotation',90,'FontSize',12);%
        end

        yyaxis right
        scatter(categorical(SA_Q_varnames(selidx),SA_Q_varnames(selidx)),SA_Q_delta_num(selidx,1),50,cl_RP*.6,'d','filled');%,'MarkerFaceAlpha',baralpha)
        hold on
        scatter(categorical(SA_Q_varnames(selidx),SA_Q_varnames(selidx)),SA_Q_delta_num(selidx,1)+SA_Q_delta_num(selidx,2),50,cl_DP*.6,'d','filled');%,'MarkerFaceAlpha',baralpha) %manually stack scatter
        ylim([-210 140])
        if jj==3 % add right ylabel to right col
            ylabel('\Delta # of plants')
        end
        ax = gca;
        ax.YAxis(2).Color = 'k';
        set(gca,'XTickLabel',[])
        grid minor
        box on
        if ii==2 && jj==2
            % Add rectangle to show default
            text(1.65,100,"Default",'FontAngle','italic')
            yy=ylim;
            rectangle('Position',[1.5 yy(1)+2 1 sum(abs(yy))-5],'LineStyle','-.') %[x y w h]
        end
        if ii==1  % add title to top row
            title(sprintf('Large DP = %s',largeDP_Qs{jj}),'FontSize',12);
        elseif ii==3  % add xlabel to bottom row
            set(gca,'XTickLabel',smallDP_Qs)
            text(1.6, -135,'Small DP','FontWeight','bold','FontSize',12);%
        end
        tt=tt+1;
        runcnt=runcnt+3;
    end
end
sgtitle("SA Qs Sust - Delta change")
%legend(tot_allscen.Properties.VariableNames{[2:3,5:6]}) %,'Location','southoutside','Orientation','horizontal')%)
lgd=legend([" ", " ", planttypes ],'NumColumns',2,'Location','southoutside','Box','off');%)
lgd.Title.String="Energy      Number         Plant type";
%set()

%% FINAL SA_SustQ: Plot bar matrix for delta SA as % of total under default scenarios - secondary yaxis not needed!
% Get delta change as % of total (RP and DP)
tot_default=sum(tot_allscen{selmixsust,2:3});
tot_default_num=sum(tot_allscen{selmixsust,5:6});

SA_Q_delta_prct=SA_Q_delta/tot_default*100;
SA_Q_delta_num_prct=SA_Q_delta_num/tot_default_num*100;

tt=1;
runcnt=1;
f=figure; %subplot(3,1,2)
for ii=1:numel(RP_Qs)
    hold all
    for jj=1:numel(largeDP_Qs)
        subtightplot (3, 3, tt, [0.025 0.08],[.08 .11],[.1 .1]); %[vgap hgap], marg_h,marg_w subplot(3,3,tt);
        selidx=runcnt:runcnt+2;
        b1=bar(categorical(SA_Q_varnames(selidx),SA_Q_varnames(selidx)),SA_Q_delta_prct(selidx,:),'stacked','EdgeAlpha',0);
        for k = 1:2
            b1(k).FaceColor = mycolors{k};
            b1(k).EdgeColor = mycolors{k};
        end
        alpha(b1,0.7)
        if jj==1 % add left ylabel to left col
            if ii==2; ylabel('% change','FontSize',12); end
            text(-01,-5.5,sprintf('RP = %s',RP_Qs{ii}),'FontWeight','bold','Rotation',90,'FontSize',12,'HorizontalAlignment','center');%
        end

        hold on
        scatter(categorical(SA_Q_varnames(selidx),SA_Q_varnames(selidx)),SA_Q_delta_num_prct(selidx,1),50,cl_RP*.6,'d','filled');%,'MarkerFaceAlpha',baralpha)
        scatter(categorical(SA_Q_varnames(selidx),SA_Q_varnames(selidx)),SA_Q_delta_num_prct(selidx,1)+SA_Q_delta_num_prct(selidx,2),50,cl_DP*.6,'d','filled');%,'MarkerFaceAlpha',baralpha) %manually stack scatter
        ylim([-25 15])
        %         if jj==3 % add right ylabel to right col
        %             if ii==2;  ylabel('% change in # of plants','FontSize',12); end
        %         end

        %ax = gca;
        %ax.YAxis(2).Color = 'k';
        set(gca,'XTickLabel',[])
        grid minor
        box on
        if ii==2 && jj==2
            % Add rectangle to show default
            yy=ylim;
            text(2,yy(2)-2,"Default",'FontAngle','italic','HorizontalAlignment','center')
            %rectangle('Position',[1.5 -44 1 34+44],'LineStyle','-.')
            rectangle('Position',[1.5 yy(1)+.5 1 sum(abs(yy))-1],'LineStyle','-.') %[x y w h]
        end
        if ii==1  % add title to top row
            title(sprintf('Large DP = %s',largeDP_Qs{jj}),'FontSize',12);
        elseif ii==3  % add xlabel to bottom row
            set(gca,'XTickLabel',smallDP_Qs)
            if jj==2; xlabel('Small DP','FontWeight','bold','FontSize',12); end
        end
        tt=tt+1;
        runcnt=runcnt+3;
    end
end
sgtitle("SA Qs Sust - Delta change in %total")
%lgd=legend(tot_allscen.Properties.VariableNames{[2:3,5:6]},'NumColumns',2) %,'Location','southoutside','Orientation','horizontal')%)
lgd=legend([" ", " ", planttypes ],'NumColumns',2,'Location','southoutside','Box','off');%)
lgd.Title.String="Energy      Number         Plant type";



%% OLD SA_SustQ: Plot bar matrix for delta SA as % of RP or DP potentials Q scenarios - Stacking of this is actually not allowed!!
SA_Q_delta_prct=SA_Q_delta./tot_allscen{selmixsust,2:3}*100;
SA_Q_delta_num_prct=SA_Q_delta_num./tot_allscen{selmixsust,5:6}*100;

tt=1;
runcnt=1;
f=figure; %subplot(3,1,2)
for ii=1:numel(RP_Qs)
    hold all
    for jj=1:numel(largeDP_Qs)
        subtightplot (3, 3, tt, [0.025 0.08],[.08 .11],[.1 .1]); %[vgap hgap], marg_h,marg_w subplot(3,3,tt);
        selidx=runcnt:runcnt+2;
        b1=bar(categorical(SA_Q_varnames(selidx),SA_Q_varnames(selidx)),SA_Q_delta_prct(selidx,:),'stacked','EdgeAlpha',0);
        for k = 1:2
            b1(k).FaceColor = mycolors{k};
            b1(k).EdgeColor = mycolors{k};
        end
        alpha(b1,0.7)
        ylim([-45 35])
        if jj==1 % add left ylabel to left col
            if ii==2; ylabel('% change','FontSize',12); end
            text(-01,-5.5,sprintf('RP = %s',RP_Qs{ii}),'FontWeight','bold','Rotation',90,'FontSize',12,'HorizontalAlignment','center');%
        end

        %yyaxis right
        hold on
        scatter(categorical(SA_Q_varnames(selidx),SA_Q_varnames(selidx)),SA_Q_delta_num_prct(selidx,1),50,cl_RP*.6,'d','filled');%,'MarkerFaceAlpha',baralpha)
        scatter(categorical(SA_Q_varnames(selidx),SA_Q_varnames(selidx)),SA_Q_delta_num_prct(selidx,1)+SA_Q_delta_num_prct(selidx,2),50,cl_DP*.6,'d','filled');%,'MarkerFaceAlpha',baralpha) %manually stack scatter
        ylim([-45 35])
        %         if jj==3 % add right ylabel to right col
        %             if ii==2;  ylabel('% change in # of plants','FontSize',12); end
        %         end

        %ax = gca;
        %ax.YAxis(2).Color = 'k';
        set(gca,'XTickLabel',[])
        grid minor
        box on
        if ii==2 && jj==2
            % Add rectangle to show default
            text(2,31,"Default",'FontAngle','italic','HorizontalAlignment','center')
            %rectangle('Position',[1.5 -44 1 34+44],'LineStyle','-.')
            yy=ylim;
            rectangle('Position',[1.5 yy(1)+1 1 sum(abs(yy))-1],'LineStyle','-.') %[x y w h]
        end
        if ii==1  % add title to top row
            title(sprintf('Large DP = %s',largeDP_Qs{jj}),'FontSize',12);
        elseif ii==3  % add xlabel to bottom row
            set(gca,'XTickLabel',smallDP_Qs)
            if jj==2; xlabel('Small DP','FontWeight','bold','FontSize',12); end
        end
        tt=tt+1;
        runcnt=runcnt+3;
    end
end
sgtitle("SA Qs Sust - Delta change in %RPDP")
%lgd=legend(tot_allscen.Properties.VariableNames{[2:3,5:6]},'NumColumns',2) %,'Location','southoutside','Orientation','horizontal')%)
lgd=legend([" ", " ", planttypes ],'NumColumns',2,'Location','southoutside','Box','off');%)
lgd.Title.String="Energy      Number         Plant type";

%% Write to excel
if saveout2xl
    varnames=tot_allscen{sel_SA_Q_idx, 1};
    cl_varnames=erase(varnames ,{'SA-','-Full-Mixed-Sust-RiskAverse'});
    SA_Q_out=table(cl_varnames,varnames, SA_Q_delta,SA_Q_delta_num, SA_Q_delta_prct, SA_Q_delta_num_prct,...
        'VariableNames',{'shortRunname','FullRunname','Delta TWh/yr RP/DP','Delta # RP/DP','Delta % TWh/yr RP/DP','Delta % # RP/DP'});


    xlsfile_summary=fullfile(rootf,'output',[oxlfname '.xlsx']);
    tmp=sprintf("Last updated on: %s",datestr(datetime('now')));

    SA_Q_out(size(SA_Q_out,1)+1,1)={tmp{:}};
    writetable(SA_Q_out,xlsfile_summary,'Sheet',"SA_Q_delta_%Tot");
end