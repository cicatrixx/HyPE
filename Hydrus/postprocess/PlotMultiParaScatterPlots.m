% Create plot of Head vs Q vs GWh vs COE
% Prepare data (takes time) or load preprocessed data for main scenarios
tic 
clc
clear vars
close all
addpath(genpath(fullfile(pwd,'..','GitConnect\Hydrus\devFiles\')), ...
    genpath(fullfile(pwd,'Hydrus')))
run('myVarNames.m')

% load('G:\SurfDrive\HPmodel\output\Figs_trial\MainScenarios4R.mat', 'pcsout',...
%     'runnames','pcsouthaz',  'runnameshaz','basindata','catchments_cl')
%  pcsout=[pcsout pcsouthaz];
% runnames=[runnames runnameshaz];
write2mat=01;
animateGIF=0;
GIF_filepath = fullfile(rootf,'output','Figs_trial','cl','DistOfPlantCharacteristics_FullOnly.gif'); % Specify the output file name

%% Load data or process and save it for future use
if ~write2mat
    load(fullfile(rootf,'output','Figs_trial','MainScenarios_BasinOut.mat'))
    disp("Loaded saved output data for Basin101")
else
    %% Loop through and compile outputs for FULL 3 policy types
    runcnt=1;
    runname_prefix='R103_Energy';
    for scen={'Full', 'Remain'}%   'Full';  % 'Remain' ;  %
        for policyname={'Large','Medium','Mixed'}
            for consttype={'Tech_Fin', 'Sust_RiskAverse'}
                % Prepare exact filename for matfile
                runnames{runcnt}=strjoin([{runname_prefix},scen(:)', policyname(:)', consttype(:)'],'_');

                % load basin output data           ofname=strcat(selrun,'\\Basin101_output_do.mat');
                pcsout_basinOut{runcnt} = load(fullfile(rootf, 'output',scen{:}, continent_in,runnames{runcnt},...
                    sprintf('Basin%d_output_do.mat', nbasin)));
                pcsout_basinOut{runcnt}.runname=runnames{runcnt};
                runcnt=runcnt+1;
            end
        end
    end
    disp("Loaded output data for Basin103 - Energy ")

    %% Loop through and compile outputs for FULL Sust Hazard Rep scenarios
    runname_prefix='R103_HazRep';
    for scen={'Full'}%   'Full';  % 'Remain' ;  %
        for policyname={'Mixed'}
            for consttype={'Sust_Costbased', 'Sust_Multi-hazard'}
                % Prepare exact filename for matfile
                runnames{runcnt}=strjoin([{runname_prefix},scen(:)', policyname(:)', consttype(:)'],'_');

                % load basin output data           ofname=strcat(selrun,'\\Basin101_output_do.mat');
                pcsout_basinOut{runcnt} = load(fullfile(rootf, 'output',scen{:}, continent_in,runnames{runcnt},...
                    sprintf('Basin%d_output_do.mat', nbasin)));
                pcsout_basinOut{runcnt}.runname=runnames{runcnt};
                runcnt=runcnt+1;
            end
        end
    end
    disp("Loaded output data for Basin103 - hazard rep")
    runnames=runnames'

    %% Combine design params for DP and RP project
    runnames_cl=erase(strrep(extractAfter(runnames,'Energy_'),'_','-'),["-Fin"])';
    myscens=1:length(runnames);

    for i=1:length(myscens)
        mydata=pcsout_basinOut{myscens(i)};
        %From model code, I know that DP data is above RP data COEAll = [horzcat(PCOEend{:})';COETotRD]; PnetAll = [horzcat(PPnetend{:})';RDPnet];
        COEalls{i}=mydata.COEAlls;
        SysIDs{i}=mydata.SysIDs;
        Pnetalls_GWhs{i}=[horzcat(mydata.PPnetend{:})'; mydata.RDPnet];  % these are in kWh in functions but converted to GWh in main script

        Qdesigns{i}=[horzcat(mydata.QDesignPinletMinend{:})'; mydata.Q_RD_design(:)];
        Heads{i}=[horzcat(mydata.HeadPminend{:})'; mydata.OptDH(:)];
        MWs{i}=[horzcat(mydata.PPMinend{:})';   mydata.RDP(:)]/1e6; % W to MW
        inValid=isnan(COEalls{i});
        Pnetalls_GWhs{i}(inValid)=nan;
        Qdesigns{i}(inValid)=nan;
        Heads{i}(inValid)=nan;
        MWs{i}(inValid)=nan;
    end
    %% Save all files
    save(fullfile(rootf,'output','Figs_trial','MainScenarios_BasinOut.mat'), ...
        'runnames', 'runnames_cl', 'COEalls', 'SysIDs', 'Pnetalls_GWhs',...
        'Qdesigns', 'Heads','MWs')

    % Can save all structure but this is not really used again so i can
    % skip for now
%     save(fullfile(rootf,'output','Figs_trial','MainScenarios_BasinOutAll.mat'), ...
%         'pcsout_basinOut','runnames', '-v7.3' )
    disp("Saved output data for Basin101 - hazard rep")

end

%% FINAL: Overlay 3 scenarios Mixed-Tech, Mixed-Costbased, Mixed-SustRiskAverse (COE)
runnames_cl=erase(strrep(extractAfter(runnames,'Full_'),'_','-'),["-Fin"])';
tmp=flipud(cbrewer2( 'Greys', 4));
mygray=tmp(1,:); '#636363';
myorange=cmap3_mainenergy(2,:); %[217,95,14]/265; %'#d95f0e';
mygreen=[brighten(cmap3_mainenergy(3,:),-0.5); brighten(cmap3_mainenergy(3,:),0.6)];

fig=figure('Position', get(0, 'Screensize'),'color','w');
% Plot mixed tech in orange
i=5;
selDPs=SysIDs{i}==1;
selRPs=SysIDs{i}==2;
hold all
s1=scatter(Qdesigns{i}(selDPs),Heads{i}(selDPs),25+Pnetalls_GWhs{i}(selDPs),myorange,'filled','MarkerFaceAlpha',0.3', 'DisplayName',planttypes{2});% ,'MarkerEdgeColor',0.7*[1 1 1]);% 'MarkerEdgeColor','k'    'MarkerFaceColor',[0 .8 .8]) ;%,'x','Color',copts2(3,:))
s2=scatter(Qdesigns{i}(selRPs),Heads{i}(selRPs),25+Pnetalls_GWhs{i}(selRPs),myorange,'^','filled','MarkerFaceAlpha',0.3','DisplayName',planttypes{1});%,'MarkerEdgeColor',0.7*[1 1 1]) ;%,'x','Color',copts2(3,:))
%title('Relation between Plant type (symbol type), GWh (symbol size), COE (symbol color), H(y axis) and Q (x axis)')
%subtitle(sprintf("SCENARIO: %s",runnames_cl{i}),'Interpreter', 'none')
%set(gca, 'YScale', 'log','XScale', 'log','ColorScale','log')
%grid on; box on

% Plot mixed sust cost in gray
i=13;
selDPs=SysIDs{i}==1;
selRPs=SysIDs{i}==2;
%hold all
s3=scatter(Qdesigns{i}(selDPs),Heads{i}(selDPs),25+Pnetalls_GWhs{i}(selDPs),mygray,'filled','MarkerFaceAlpha',0.3', 'DisplayName',planttypes{2});% ,'MarkerEdgeColor',0.7*[1 1 1]'MarkerEdgeColor','k'    'MarkerFaceColor',[0 .8 .8]) ;%,'x','Color',copts2(3,:))
s4=scatter(Qdesigns{i}(selRPs),Heads{i}(selRPs),25+Pnetalls_GWhs{i}(selRPs),mygray,'^','filled','MarkerFaceAlpha',0.3','DisplayName',planttypes{1}) ;% ,'MarkerEdgeColor',0.7*[1 1 1],'x','Color',copts2(3,:))
% title('Relation between Plant type (symbol type), GWh (symbol size), COE (symbol color), H(y axis) and Q (x axis)')
% subtitle(sprintf("SCENARIO: %s",runnames_cl{i}),'Interpreter', 'none')
% set(gca, 'YScale', 'log','XScale', 'log','ColorScale','log')
% grid on; box on

% Plot mixed sust hazardbased w blackoutline = this has strong overlap w
% the cost based case, only few big triangles are left out
% i=14;
% selDPs=SysIDs{i}==1;
% selRPs=SysIDs{i}==2;
% hold all
% s31=scatter(Qdesigns{i}(selDPs),Heads{i}(selDPs),25+Pnetalls_GWhs{i}(selDPs),'k','MarkerFaceAlpha',0,'LineWidth',1.5, 'DisplayName',planttypes{2});% ,'MarkerEdgeColor',0.7*[1 1 1] 'MarkerEdgeColor','k'    'MarkerFaceColor',[0 .8 .8]) ;%,'x','Color',copts2(3,:))
% s32=scatter(Qdesigns{i}(selRPs),Heads{i}(selRPs),25+Pnetalls_GWhs{i}(selRPs),'k','^','MarkerFaceAlpha',0,'LineWidth',1.5,'DisplayName',planttypes{1}) ;%,'MarkerEdgeColor',0.7*[1 1 1] ,'x','Color',copts2(3,:))

% Plot mixed sust riskaverse w greenoutline
i=6;
selDPs=SysIDs{i}==1;
selRPs=SysIDs{i}==2;
hold all
%s5=scatter(Qdesigns{i}(selDPs),Heads{i}(selDPs),25+Pnetalls_GWhs{i}(selDPs),COEalls{i}(selDPs),'filled','MarkerFaceAlpha',0.6', 'DisplayName',planttypes{2});% ,'MarkerEdgeColor',0.7*[1 1 1] 'MarkerEdgeColor','k'    'MarkerFaceColor',[0 .8 .8]) ;%,'x','Color',copts2(3,:))
%s6=scatter(Qdesigns{i}(selRPs),Heads{i}(selRPs),25+Pnetalls_GWhs{i}(selRPs),COEalls{i}(selRPs),'^','filled','MarkerFaceAlpha',0.6','DisplayName',planttypes{1}) ;%,'MarkerEdgeColor',0.7*[1 1 1] ,'x','Color',copts2(3,:))
s5=scatter(Qdesigns{i}(selDPs),Heads{i}(selDPs),25+Pnetalls_GWhs{i}(selDPs),mygreen(1,:),'MarkerFaceAlpha',0,'LineWidth',1, 'DisplayName',planttypes{2});% ,'MarkerEdgeColor',0.7*[1 1 1] 'MarkerEdgeColor','k'    'MarkerFaceColor',[0 .8 .8]) ;%,'x','Color',copts2(3,:))
s6=scatter(Qdesigns{i}(selRPs),Heads{i}(selRPs),25+Pnetalls_GWhs{i}(selRPs),mygreen(1,:),'^','MarkerFaceAlpha',0,'LineWidth',1,'DisplayName',planttypes{1}) ;%,'MarkerEdgeColor',0.7*[1 1 1] ,'x','Color',copts2(3,:))

axis square; lgdbox=[0.370486121650577,0.114876718379524,0.357812489460533,0.084029225317571];
set(gca, 'YScale', 'log','XScale', 'log','ColorScale','log',...
    'FontName', 'Segoe UI','FontSize',14,'Layer','top');
set(gca,'yticklabels',num2str(yticks'),'xticklabels',num2str(xticks'))
ylabel('Head in m','FontWeight','bold')
xlabel('Design discharge in m^3/s','FontWeight','bold')
grid on; box on

%mycbar("Cost of production (USD/kWh)")
title('Relation between Plant type (symbol type), GWh (symbol size), COE (symbol color), H(y axis) and Q (x axis)')
subtitle(sprintf("SCENARIO: %s",runnames_cl{i}),'Interpreter', 'none')
%colormap(cbrewer2('YlGnBu')) ;colormap(viridis)

%lgdbox=[0.547569454983911,0.100904662491301,0.357812489460533,0.084029225317571];
lgd=legend([s2 s1 s4 s3 s6 s5],[repmat("                            .",1,4) strcat("                    .",planttypes)],...
    'NumColumns',3,'Location','southeast','EdgeColor','none','Position',lgdbox);
title(lgd,"Technical    Sustainable: Cost based    Sustainable: Risk averse          PLANT TYPE          ")
annotation('rectangle', lgdbox, 'LineWidth',0.3);

%% Individual: ENERGY SCENS 5-in-one H vs Q symbol size= GWh, color =COE for RP and DP type
% Combine design params for DP and RP project
i=1 ; % Full large tech

figure;
%subplot(1,2,1);
selDPs=SysIDs{i}==1;
selRPs=SysIDs{i}==2;
hold all
scatter(Qdesigns{i}(selDPs),Heads{i}(selDPs),25+Pnetalls_GWhs{i}(selDPs),COEalls{i}(selDPs),'filled','MarkerFaceAlpha',0.6', 'DisplayName',planttypes{2},'MarkerEdgeColor',0.7*[1 1 1]);% 'MarkerEdgeColor','k'    'MarkerFaceColor',[0 .8 .8]) ;%,'x','Color',copts2(3,:))
scatter(Qdesigns{i}(selRPs),Heads{i}(selRPs),25+Pnetalls_GWhs{i}(selRPs),COEalls{i}(selRPs),'^','filled','MarkerFaceAlpha',0.6','DisplayName',planttypes{1},'MarkerEdgeColor',0.7*[1 1 1]) ;%,'x','Color',copts2(3,:))
set(gca, 'YScale', 'log','XScale', 'log','ColorScale','log')
ylabel('Head in m')
xlabel('Design discharge in m^3/s')
grid on; box on
legend()
mycbar("Cost of production (USD/kWh)")
title('Relation between Plant type (symbol type), GWh (symbol size), COE (symbol color), H(y axis) and Q (x axis)')
subtitle(sprintf("SCENARIO: %s",runnames_cl{myscens(i)}),'Interpreter', 'none')
colormap(viridis)

%% GOOD: Plot GIF of distribution for 6 energy scenarios
myscens=[1:6 13 14 6];
nplots=length(myscens);
runnames_cl=strrep(erase(runnames,["R101_","Energy_","HazRep_","_Fin","Energy_"]),'_','-')';

% Create Plots
%Set up fix for making GIF
fig=figure('Position', [.5,.5 1000 1000]);
axes('Units', 'normalized', 'Position', [0 0 1 1])

for m=1:nplots
    %Plot fig
    i=myscens(m);
    selDPs=SysIDs{i}==1;
    selRPs=SysIDs{i}==2;
    clf % use disable clf to ID axis limit first
    hold all
    s1=scatter(Qdesigns{i}(selDPs),Heads{i}(selDPs),25+Pnetalls_GWhs{i}(selDPs),COEalls{i}(selDPs),'filled','MarkerFaceAlpha',0.6', 'DisplayName',planttypes{2});% ,'MarkerEdgeColor',0.7*[1 1 1] 'MarkerEdgeColor','k'    'MarkerFaceColor',[0 .8 .8]) ;%,'x','Color',copts2(3,:))
    s2=scatter(Qdesigns{i}(selRPs),Heads{i}(selRPs),25+Pnetalls_GWhs{i}(selRPs),COEalls{i}(selRPs),'^','filled','MarkerFaceAlpha',0.6','DisplayName',planttypes{1}) ;%,'MarkerEdgeColor',0.7*[1 1 1] ,'x','Color',copts2(3,:))
    
    %Labels
    ylabel('Head in m','FontWeight','bold')
    xlabel('Design discharge in m^3/s','FontWeight','bold')
    grid on; box on
    legend([s2 s1],planttypes,'Position',[0.626704633527781,0.166333334763845,0.21278720796406,0.052499998569488])
    mycbar("Cost of production (USD/kWh)")
    title(["Relation between plant type (symbol type), GWh (symbol size),"; "COE (symbol color), head (y axis) and design discharge (x axis)"],'Position',[26 2300])
    subtitle(sprintf("SCENARIO: %s",runnames_cl{i}),'Interpreter', 'none','HorizontalAlignment','center')%,'Position',[0.1 1600])
    colormap(viridis)

    % set limits for plot
    axis square
    caxis([0.03 320])
    xlim([.1 5000])
    ylim([1 1500])
    %change axis labels

    xticks([.1 , 1, 10 ,100, 1000])
    xticklabels([.1 , 1, 10 ,100, 1000])
    yticks([1, 10 ,100, 1000])
    yticklabels([1, 10 ,100, 1000])

set(gca, 'YScale', 'log','XScale', 'log','ColorScale','log',...
        'FontName', 'Segoe UI','FontSize',14,'Layer','top');
    %Save fig to frame
    drawnow
    frame = getframe(fig);
    im{m} = frame2im(frame);
end

%     %% Show plots as subplots
%     figure;
%     for m=1:nplots
%         subplot(ceil(nplots/2),2,m)
%         %figure
%         imshow(im{m});
%     end

%% Animate Plots
if animateGIF
    delaytsec=2;
    for selstations=1:nplots
        [A,map] = rgb2ind(im{selstations},256);
        if selstations == 1
            imwrite(A,map,GIF_filepath,'gif','LoopCount',Inf,'DelayTime',delaytsec);
        else
            imwrite(A,map,GIF_filepath,'gif','WriteMode','append','DelayTime',delaytsec);
        end
    end
end

%orient('landscape');print('tmp.pdf','-dpdf','-bestfit')

