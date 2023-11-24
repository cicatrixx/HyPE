% Prepare data for cost plots
% Add one plot w 3x2 energy + 2 hazard rep cases
% Add plot w sub-basin wise distinction of mixed case

% Load data for main scenarios - same as gettotals... then get sub-basin wise sums
clc
clearvars
close all
run('myVarNames.m')

load('G:\SurfDrive\HPmodel\output\Figs_trial\MainScenarios.mat', 'pcsout',...
    'runnames','pcsouthaz',  'runnameshaz','basindata','catchments_cl')
pcsout=[pcsout pcsouthaz];
runnames=[runnames runnameshaz];

if finalfig
    %% GOOD: Loop FULL + Hazard Cost curve w points for DP and RP - assumes tech and sust scenarios are back to back
    fig=figure('Position', get(0, 'Screensize'),'color','w');
    
    subplot(1,3,1:2)
    runnames_cl=strrep(strrep(extractAfter(runnames,'Full_'),'_Tech_Fin','-Technical'),'_Sust_','-Sustainable-'); %,["Large-", "Medium-", "Mixed-"]);
    myscens=[2 1 4 3 6 5 13 14];
    copts=[];
    %create my color palette with one shade lighter version for Tech case
    for i=1:3 % 3 scens % create lighter shade for each scenario
        copts=[copts; cmap3_mainenergy(i,:);  brighten(cmap3_mainenergy(i,:),0.6);  ];
    end

    % Add greys for hazard
    copts=[copts; flipud(cbrewer2( 'Greys', 4))];
    c=1;
    mylegend={};
    mylegend_scname={};
    for scn= myscens
        scol=copts(c,:);c=c+1; %bcoz i may not go serially
        plotCostCurve(pcsout{scn}.COEAlls, pcsout{scn}.PnetAlls, pcsout{scn}.SysIDs,scol,'-',.8,runnames_cl{scn});
        %    plotCostCurveScatter(pcsout{i}.COEAlls, pcsout{i}.PnetAlls, pcsout{i}.SysIDs,scol,1);
        %sgtitle(selrun,'Interpreter', 'none' )
        mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
        mylegend_scname= [mylegend_scname(:)' {''} {''} runnames_cl(scn) ];
                end
    % Add lines for potential cutoffs
    %yline(Techlim,'LineStyle','--','LineWidth',1.25); %,'DisplayName','Technical potential limit'); %Tech potential
    %text(10,Techlim*1.15, sprintf("\\it Technical potential\n<= %0.2f USD/kWh",Techlim),'FontSize',12)
    yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(240,0.05, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)
    
    xlim([0,305])
    ylabel('Per unit production cost (USD2010 per KWh)',FontWeight='bold')
    xlabel('Cumulative Energy (TWh)',FontWeight='bold')
    set(gca,'YScale','log')
    set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
    %set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
    grid on; box on
    %lgd=legend(mylegend,'NumColumns',2,'Orientation','horizontal','Location','best');
    %lgd=legend(mylegend_scname,'NumColumns',3,'Orientation','vertical', 'Interpreter','none','Location','north'); %
    lgdbox=[0.301 0.844 0.32239582768331 0.071503129657475];
    lgd=legend([planttypes {'   '} repmat([{''} {''} {'   '}],1,3) {''} {''} "   Sustainable" {''} {''} "   Technical" {''} {''} "   Cost-Based" {''} {''} "   Multi-hazard" ],...
        'NumColumns',5,'Orientation','vertical','Location','northoutside', 'Interpreter','none',...
        'EdgeColor','none','Position',lgdbox);
    annotation('rectangle',...
    lgdbox,...
    'LineWidth',0.3);

    title(lgd,sprintf('         PLANT TYPE            %s %s %s    SCENARIO            RISK TYPE',searchtypes{:}))
    % savefig(fullfile(outfig_fldr,'CostCurve_AllMajor.fig'))
    % savefig(fullfile(outfig_fldr,'CostCurve2_ForEGU2021.fig'))
    % export_fig(fullfile(outfig_fldr,'CostCurve2_ForEGU2021.jpg'))
    myylim=ylim;

    % GOOD: Sub-basinwise cost curves for FULL Mixed sust
    %figure
    subplot(1,3,3)
    cmap8=cmap8_wong;
    scn=6; % the Full mixed sust scenario
    addfill=0;
    %subplot(4,1,2:4)
    mylegend={};
    mylegend_basname={};
    c=1;
    for basID= 101:108
        scol=cmap8(c,:);c=c+1; %bcoz i may not go serially
        selbas=pcsout{scn}.subbasinIDs==basID;
        plotCostCurveScatter(pcsout{scn}.COEAlls(selbas), pcsout{scn}.PnetAlls(selbas), pcsout{scn}.SysIDs(selbas),scol,addfill,'-',basindata.basinnames{basID-100});
        mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
        mylegend_basname= [mylegend_basname(:)' {''} {''} [basindata.basinnames(basID-100) ]];
    end

    % Add lines for potential cutoffs
    %yline(Techlim,'LineStyle','--','LineWidth',1.25); %,'DisplayName','Technical potential limit'); %Tech potential
    %text(10,Techlim*1.15, sprintf("\\it Technical potential\n<= %0.2f USD/kWh",Techlim),'FontSize',12)
    yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(24,0.05, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)
    xlabel('Cumulative Energy (TWh)',FontWeight='bold')
    ylabel('Per unit production cost (USD2010 per KWh)',FontWeight='bold')

    ylim(myylim)
    set(gca,'YScale','log')
    set(gca,'yticklabels',num2str(yticks')) %use number labels w full numbers and not the exponential form!
    grid on; box on
    sgtitle("Mixed sust")

    % subcatchment map overlaid
    %subplot(4,1,1)
%    axes('Position', [1-.3-.095, .75, .2, .2]) % for independent fig
    axes('Position', [1-.277, .8, .15, .15]) % for subplot case
    h=imagescnan(catchments_cl);
    colormap(cmap8)
    set(h, 'AlphaData', myalpha*(catchments_cl>0))
    
    % add names for subbasins
    text(bas_x(1:8),bas_y(1:8),basindata.basinnames,"HorizontalAlignment","center","FontAngle","italic")

%     cb=colorbar('Ticks', basindata.basinIDs,'TickLabels',basindata.basinnames,'Direction','reverse','Location','westoutside');
%     % Make the colorbar transparent - Need to do this manually
%     cdata = cb.Face.Texture.CData;
%     cdata(end,:) = uint8(myalpha * cdata(end,:));
%     cb.Face.Texture.ColorType = 'truecoloralpha';
%     cb.Face.Texture.CData = cdata;
else
    %% GOOD: Loop FULL and REMAIN Energy Only Cost curve w points for DP and RP - assumes tech and sust scenarios are back to back
    runnames_cl=erase(strrep(extractAfter(runnames,'Energy_'),'_','-'),["Large-", "Medium-", "Mixed-"]);
    nplot=length(runnames);
    myscens=[2     1     8     7 4     3    10     9    6     5    12    11]; %reorder how scenarios are run

    % selected palettes for color blind friendly
    %copts=cbrewer2( 'Paired', 6); %linspecer(numel(myscens)/2,'qualitative');%['#404040'; '#f4a582';'#ca0020' ]; %parula(3);
    copts=[flipud(cbrewer2( 'Purples', 4));flipud(cbrewer2( 'Oranges', 4));flipud(cbrewer2( 'Greens', 4))];
    copts=brighten(copts,.2); %lighten the shades

    c=1;
    figure
    mylegend={};
    mylegend_scname={};
    for scn= myscens
        scol=copts(c,:);c=c+1; %bcoz i may not go serially
        plotCostCurve(pcsout{scn}.COEAlls, pcsout{scn}.PnetAlls, pcsout{scn}.SysIDs,scol);
        %sgtitle(selrun,'Interpreter', 'none' )
        mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
        mylegend_scname= [mylegend_scname(:)' {''} {''} [runnames_cl(scn) ]];
    end

    % Add lines for potential cutoffs
    %yline(Techlim,'LineStyle','--','LineWidth',1.25); %,'DisplayName','Technical potential limit'); %Tech potential
    %text(10,Techlim*1.15, sprintf("\\it Technical potential\n<= %0.2f USD/kWh",Techlim),'FontSize',12)
    yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(300,costlim*1.15, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)
    xlabel('Cumulative Energy (TWh)')
    ylabel('Per unit production cost (USD2010 per KWh)')

    % Format yaxis to prevent scientific notation : https://nl.mathworks.com/matlabcentral/answers/367567-how-can-i-control-the-format-of-the-tick-label-numbers-when-using-a-log-yscale-with-a-large-range
    %ylim([0 1000])
    % yticks(0:100:500)
    % yticklabels(num2str((0.01:100:500)'))
    %ytickformat('%,.0g')
    %ax = gca;
    %ax.YAxis.Exponent = 0;
    %ax.YScale = 'log';
    set(gca,'YScale','log')
    set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
    grid on; box on
    %lgd=legend(mylegend,'NumColumns',3,'Orientation','horizontal','Location','best');
    %lgd=legend(mylegend_scname,'NumColumns',3,'Orientation','vertical','Location','best', 'Interpreter','none');

    %manual legend
    lgd=legend([repmat([{''} {''} {'   '}],1,8) {''} {''} "  Full - Sustainable" {''} {''} "  Full - Technical" {''} {''} "  Remaining - Sustainable" {''} {''} "  Remaining - Technical" ],...
        'NumColumns',3,'Orientation','vertical','Location','northoutside', 'Interpreter','none',...
        'EdgeColor','none');
    title(lgd,'Large    Medium    Mixed           SCENARIO          ')
    %
    % savefig(fullfile(outfig_fldr,'CostCurve_AllMajor.fig'))
    % savefig(fullfile(outfig_fldr,'CostCurve2_ForEGU2021.fig'))
    % export_fig(fullfile(outfig_fldr,'CostCurve2_ForEGU2021.jpg'))

    %% GOOD: Loop FULL Cost curve w points for DP and RP - assumes tech and sust scenarios are back to back
    myscens=[2 1 4 3 6 5];
    copts=[flipud(cbrewer2( 'Purples', 2));flipud(cbrewer2( 'Oranges', 2));flipud(cbrewer2( 'Greens', 2))];
    copts=brighten(copts,.2);
    c=1;
    figure
    mylegend={};
    mylegend_scname={};
    for scn= myscens
        scol=copts(c,:);c=c+1; %bcoz i may not go serially
        plotCostCurve(pcsout{scn}.COEAlls, pcsout{scn}.PnetAlls, pcsout{scn}.SysIDs,scol);
        %sgtitle(selrun,'Interpreter', 'none' )
        mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
        mylegend_scname= [mylegend_scname(:)' {''} {''} runnames_cl(scn) ];
    end

    % Add lines for potential cutoffs
    %yline(Techlim,'LineStyle','--','LineWidth',1.25); %,'DisplayName','Technical potential limit'); %Tech potential
    %text(10,Techlim*1.15, sprintf("\\it Technical potential\n<= %0.2f USD/kWh",Techlim),'FontSize',12)
    yline(costlim,'LineStyle','--','Color','k','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(300,costlim*1.15, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)

    ylabel('Per unit production cost (USD2010 per KWh)')
    xlabel('Cumulative Energy (TWh)')
    set(gca,'YScale','log')
    set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
    grid on; box on
    %lgd=legend(mylegend,'NumColumns',2,'Orientation','horizontal','Location','best');
    lgd=legend(mylegend_scname,'NumColumns',3,'Orientation','vertical','Location','best', 'Interpreter','none');
    lgd=legend([planttypes {'   '} repmat([{''} {''} {'   '}],1,3) {''} {''} "   Sustainable" {''} {''} "   Technical"],...
        'NumColumns',4,'Orientation','vertical','Location','northoutside', 'Interpreter','none',...
        'EdgeColor','none');
    title(lgd,'         PLANT TYPE           Large   Medium   Mixed   SCENARIO')
    % savefig(fullfile(outfig_fldr,'CostCurve_AllMajor.fig'))
    % savefig(fullfile(outfig_fldr,'CostCurve2_ForEGU2021.fig'))
    % export_fig(fullfile(outfig_fldr,'CostCurve2_ForEGU2021.jpg'))


    %% GOOD: Loop FULL + Hazard Cost curve w points for DP and RP - assumes tech and sust scenarios are back to back
    runnames_cl=strrep(extractAfter(runnames,'Full_'),'_','-'); %,["Large-", "Medium-", "Mixed-"]);
    myscens=[2 1 4 3 6 5 13 14];
    %create my color palette
    copts=[flipud(cbrewer2( 'Purples', 2));flipud(cbrewer2( 'Oranges', 2));flipud(cbrewer2( 'Greens', 2))];
    copts=[copts; flipud(cbrewer2( 'Greys', 4))];

    c=1;
    figure
    mylegend={};
    mylegend_scname={};
    for scn= myscens
        scol=copts(c,:);c=c+1; %bcoz i may not go serially
        plotCostCurve(pcsout{scn}.COEAlls, pcsout{scn}.PnetAlls, pcsout{scn}.SysIDs,scol);
        %    plotCostCurveScatter(pcsout{i}.COEAlls, pcsout{i}.PnetAlls, pcsout{i}.SysIDs,scol,1);

        %sgtitle(selrun,'Interpreter', 'none' )
        mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
        mylegend_scname= [mylegend_scname(:)' {''} {''} runnames_cl(scn) ];
    end

    % Add lines for potential cutoffs
    %yline(Techlim,'LineStyle','--','LineWidth',1.25); %,'DisplayName','Technical potential limit'); %Tech potential
    %text(10,Techlim*1.15, sprintf("\\it Technical potential\n<= %0.2f USD/kWh",Techlim),'FontSize',12)
    yline(costlim,'LineStyle','--','Color','k','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(300,costlim*1.15, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)

    ylabel('Per unit production cost (USD2010 per KWh)')
    xlabel('Cumulative Energy (TWh)')
    set(gca,'YScale','log')
    set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
    grid on; box on
    %lgd=legend(mylegend,'NumColumns',2,'Orientation','horizontal','Location','best');
    lgd=legend(mylegend_scname,'NumColumns',3,'Orientation','vertical','Location','best', 'Interpreter','none');
    lgd=legend([planttypes {'   '} repmat([{''} {''} {'   '}],1,3) {''} {''} "   Sustainable" {''} {''} "   Technical" {''} {''} "   Cost-Based" {''} {''} "   Multi-hazard" ],...
        'NumColumns',5,'Orientation','vertical','Location','northoutside', 'Interpreter','none',...
        'EdgeColor','none');
    title(lgd,'         PLANT TYPE            Large    Medium   Small    SCENARIO            RISK TYPE')
    % savefig(fullfile(outfig_fldr,'CostCurve_AllMajor.fig'))
    % savefig(fullfile(outfig_fldr,'CostCurve2_ForEGU2021.fig'))
    % export_fig(fullfile(outfig_fldr,'CostCurve2_ForEGU2021.jpg'))


    %% GOOD: Loop FULL+REMAIN+Hazard Cost curve w points for DP and RP - assumes tech and sust scenarios are back to back
    runnames_cl=strrep(extractAfter(runnames,'R101_'),'_','-');
    nplot=length(runnames);
    myscens=[2     1     8     7 4     3    10     9    6     5    12    11 13 14]; %reorder how scenarios are run

    % Create color scheme
    cmapTFS=[flipud(cbrewer2( 'Purples', 4));flipud(cbrewer2( 'Oranges', 4));flipud(cbrewer2( 'Greens', 4))];

    % selected palettes for color blind friendly
    %copts=cbrewer2( 'Paired', 6); %linspecer(numel(myscens)/2,'qualitative');%['#404040'; '#f4a582';'#ca0020' ]; %parula(3);
    copts=[flipud(cbrewer2( 'Purples', 4));flipud(cbrewer2( 'Oranges', 4));flipud(cbrewer2( 'Greens', 4))];
    copts=[copts; flipud(cbrewer2( 'Greys', 4))];
    %copts=brighten(copts,.2); %lighten the shades
    % Reorder
    %myscens=[8 7 10 9 12 11 2 1  4 3 6 5 13 14]; %reorder how scenarios are run
    %copts=[copts([3 4 7 8 11 12 1 2 5 6 9 10] ,:); flipud(cbrewer2( 'Greys', 4))];

    c=1;
    figure
    mylegend={};
    mylegend_scname={};
    for scn= myscens
        scol=copts(c,:);c=c+1; %bcoz i may not go serially
        plotCostCurve(pcsout{scn}.COEAlls, pcsout{scn}.PnetAlls, pcsout{scn}.SysIDs,scol);
        %sgtitle(selrun,'Interpreter', 'none' )
        mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
        mylegend_scname= [mylegend_scname(:)' {''} {''} [runnames_cl(scn) ]];
    end

    % Add lines for potential cutoffs
    %yline(Techlim,'LineStyle','--','LineWidth',1.25); %,'DisplayName','Technical potential limit'); %Tech potential
    %text(10,Techlim*1.15, sprintf("\\it Technical potential\n<= %0.2f USD/kWh",Techlim),'FontSize',12)
    yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(300,costlim*1.15, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)
    xlabel('Cumulative Energy (TWh)')
    ylabel('Per unit production cost (USD2010 per KWh)')

    % Format yaxis to prevent scientific notation : https://nl.mathworks.com/matlabcentral/answers/367567-how-can-i-control-the-format-of-the-tick-label-numbers-when-using-a-log-yscale-with-a-large-range
    %ylim([0 1000])
    % yticks(0:100:500)
    % yticklabels(num2str((0.01:100:500)'))
    %ytickformat('%,.0g')
    %ax = gca;
    %ax.YAxis.Exponent = 0;
    %ax.YScale = 'log';
    set(gca,'YScale','log')
    set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
    grid on; box on
    %lgd=legend(mylegend,'NumColumns',3,'Orientation','horizontal','Location','best');
    lgd=legend(mylegend_scname,'NumColumns',4,'Orientation','vertical','Location','best', 'Interpreter','none');
    lgd=legend([repmat([{''} {''} {'   '}],1,8) {''} {''} "  Full - Sustainable      " {''} {''} "  Full - Technical      " {''} {''} "  Remaining - Sustainable  " {''} {''} "  Remaining - Technical" {''} {''} "   Cost-Based" {''} {''} "   Multi-hazard" ],...
        'NumColumns',4,'Orientation','vertical','Location','northoutside', 'Interpreter','none',...
        'EdgeColor','none');
    title(lgd,'Large    Medium    Mixed         SCENARIO                    RISK TYPE    ')

    % savefig(fullfile(outfig_fldr,'CostCurve_AllMajor.fig'))
    % savefig(fullfile(outfig_fldr,'CostCurve2_ForEGU2021.fig'))
    % export_fig(fullfile(outfig_fldr,'CostCurve2_ForEGU2021.jpg'))

    %% GOOD: Cost curves w changing symbol sizes FULL+REMAIN+Hazard
    myscens=[2     1     8     7 4     3    10     9    6     5    12    11 13 14]; %reorder how scenarios are run
    addfill=0;
    figure
    mylegend={};
    mylegend_scname={};
    c=1;

    for scn= myscens
        scol=copts(c,:);c=c+1; %bcoz i may not go serially
        plotCostCurveScatter(pcsout{scn}.COEAlls, pcsout{scn}.PnetAlls, pcsout{scn}.SysIDs,scol,addfill);
        %sgtitle(selrun,'Interpreter', 'none' )
        mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
        mylegend_scname= [mylegend_scname(:)' {''} {''} [runnames_cl(scn) ]];
    end

    % Add lines for potential cutoffs
    %yline(Techlim,'LineStyle','--','LineWidth',1.25); %,'DisplayName','Technical potential limit'); %Tech potential
    %text(10,Techlim*1.15, sprintf("\\it Technical potential\n<= %0.2f USD/kWh",Techlim),'FontSize',12)
    yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(300,costlim*1.15, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)
    xlabel('Cumulative Energy (TWh)')
    ylabel('Per unit production cost (USD2010 per KWh)')

    %ax = gca;
    %ax.YAxis.Exponent = 0;
    %ax.YScale = 'log';
    set(gca,'YScale','log')
    set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
    grid on; box on
    lgd=legend([repmat([{''} {''} {'   '}],1,8) {''} {''} "  Full - Sustainable      " {''} {''} "  Full - Technical      " {''} {''} "  Remaining - Sustainable  " {''} {''} "  Remaining - Technical" {''} {''} "   Cost-Based" {''} {''} "   Multi-hazard" ],...
        'NumColumns',4,'Orientation','vertical','Location','northoutside', 'Interpreter','none',...
        'EdgeColor','none');
    title(lgd,'Large    Medium    Mixed         SCENARIO                    RISK TYPE    ')

    %% Sub-basinwise cost curves for FULL Mixed tech-econ
    %cmap8=inferno(10);
    %cmap8=cmap8(2:9,:);

    scn=5; % the Full mixed tech-econ scenario
    addfill=0;
    figure
    subplot(4,1,2:4)

    mylegend={};
    mylegend_basname={};
    c=1;
    for basID= 101:108
        scol=cmap8(c,:);c=c+1; %bcoz i may not go serially
        selbas=pcsout{scn}.subbasinIDs==basID;
        plotCostCurveScatter(pcsout{scn}.COEAlls(selbas), pcsout{scn}.PnetAlls(selbas), pcsout{scn}.SysIDs(selbas),scol,addfill);
        mylegend= [mylegend(:)' {''} {''} {char(63+c)} ];
        mylegend_basname= [mylegend_basname(:)' {''} {''} [basindata.basinnames(basID-100) ]];
    end

    % Add lines for potential cutoffs
    %yline(Techlim,'LineStyle','--','LineWidth',1.25); %,'DisplayName','Technical potential limit'); %Tech potential
    %text(10,Techlim*1.15, sprintf("\\it Technical potential\n<= %0.2f USD/kWh",Techlim),'FontSize',12)
    yline(costlim,'LineStyle','--','Color','b','LineWidth',1.25); %,'DisplayName','Economic potential limit'); %Economic potential
    text(300,costlim*1.15, sprintf("\\it Finanical potential\n<= %0.2f USD/kWh",costlim),'FontSize',12)
    xlabel('Cumulative Energy (TWh)')
    ylabel('Per unit production cost (USD2010 per KWh)')

    %ax = gca;
    %ax.YAxis.Exponent = 0;
    %ax.YScale = 'log';
    set(gca,'YScale','log')
    set(gca,'FontName','franklin gothic medium','FontSize',14,'Layer','top');
    grid on; box on
    ylim([0 80])
    %lgd=legend(mylegend_basname,'Orientation','vertical','Location','best', 'Interpreter','none');

    % lgd=legend([repmat([{''} {''} {'   '}],1,8) {''} {''} "  Full - Sustainable      " {''} {''} "  Full - Technical      " {''} {''} "  Remaining - Sustainable  " {''} {''} "  Remaining - Technical" {''} {''} "   Cost-Based" {''} {''} "   Multi-hazard" ],...
    %     'NumColumns',4,'Orientation','vertical','Location','northoutside', 'Interpreter','none',...
    %     'EdgeColor','none');
    %title(lgd,'Large    Medium    Mixed         SCENARIO                    RISK TYPE    ')


    % subcatchment map

    subplot(4,1,1)
    imagescnan(catchments_cl)
    colormap(cmap8)
    colorbar('Ticks', basindata.basinIDs,'TickLabels',basindata.basinnames,'Direction','reverse','Location','eastoutside')
    % savefig(fullfile(sprintf("%s/%s/TheoreticalPotential_SubbasinDefintion.fig",root,res)))

end