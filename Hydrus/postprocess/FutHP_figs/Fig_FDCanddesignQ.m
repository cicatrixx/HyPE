% Fig_FDCanddesignQ
% Created By    : Sanita Dhaubanjar
% Created For	: SustainIndus project WP2 Paper publication
%=========================
% Code for loading pre-processed output files from Hydrus model to create
% final plots used for publication
% Plot sub-basin wise FDC, total and % change

clc
close all
clear all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')
matout=fullfile(rootoffut, "RQatoutlets_m3s_02.mat");
load(matout, 'myFDC','myFDC_changeprct','designQ_selchanges_subbasin')
disp("Loaded files")
%% DesignQ thresholds
designQnames=["Q30", "Q40", "Q80"];
designQ_exceedances=[30, 40, 80];
designQ_selmonths=round(designQ_exceedances/100*12);

sspidx=[2:4:26];  % start idx for each ssp
sspselects=[1:3; 4:6];
monexceedances_xaxis=(1:12)/12*100;
nsubbas=size(idata_subbasin,2);
%% FINAL: Plot FDC sim hist vs fut w BANDS
subplottight = @(m,n,p) subtightplot (m, n, p, [0.05 0.05],[.11 .1],[.09 .02]); %[vgap hgap], marg_h -[low up],marg_w -[l r]

idata_subbasin=myFDC;
% For each subbasin, loop through each of the 3SSPs and plot envelope of 4 corners
for seltf=2 %
    figure(seltf*100);clf;
    sspselect=sspselects(seltf,:);

    for selsubbas=1:size(idata_subbasin,2)
        subplottight(2,4,selsubbas)
        hold all
        ii=1;
        % Plot future simulations plot
        for ssp=sspselect
            l3(ii)=plotEnvelopeMinMaxMean(squeeze(idata_subbasin(:,selsubbas,(sspidx(ssp):sspidx(ssp+1)-1))),c_ssp_main(ssp,:),3,monexceedances_xaxis);
            %scatter(1:12,squeeze(myFDC(:,selsubbas,(sspidx(ssp):sspidx(ssp+1)-1))),5,c_ssp_main(ssp,:),'o');
            ii=ii+1;
        end
        % Plot historical simulation
        hold on
        l2=plot(monexceedances_xaxis,idata_subbasin(:,selsubbas,1),"-ok");%,"DisplayName","Sim: Historical");%,'Color',[0 0 0]+.5)
        applymyplotformat( strcat(charlbl(selsubbas),basinlabels.basinnames{selsubbas}))

        % Apply labels and formatting to graph
        xline(designQ_exceedances,'LineStyle',':','Color',mygraylines) % main design discharges
        %
        if selsubbas==1
            % add designQ labels
            topy=ylim;
            text(designQ_exceedances,floor(topy(2))*[1 1 1],compose("Q%d",designQ_exceedances),'HorizontalAlignment','center','FontWeight','bold')
            % add ylabels only in left plots
            ylabel("Discharge (m^3s^{-1})")
        end

        % add xlabels only in bottom plots
        xticks([0:25:100])

        if selsubbas<=4
            xticklabels([])
        end
        if selsubbas==6
            xlabel('Flow exceedance (%)')
        end
        %decided that log scale didnt add more value
        %set(gca,'YScale','log')
        %ylim([10^1 10^4])
    end
    legend([l2 l3],modelnames(2:end),'Orientation',"horizontal")
    sgtitle(sprintf("Historical and future FDC for %s",tf_full{seltf}))
end

%% FINAL: Plot FDC sim hist vs fut w BANDS
subplottight = @(m,n,p) subtightplot (m, n, p, [0.05 0.02],[.15 .1],[.09 .02]); %[vgap hgap], marg_h -[low up],marg_w -[l r]

% For each subbasin, loop through each of the 3SSPs and plot envelope of 4 corners
idata_subbasin=myFDC_changeprct;

for seltf=2 %
    figure(seltf*10);clf;
    sspselect=sspselects(seltf,:);

    for selsubbas=1:nsubbas
        subplottight(2,4,selsubbas)
        hold all
        ii=1;
        % Plot future simulations plot
        for ssp=sspselect
            l3(ii)=plotEnvelopeMinMaxMean(squeeze(idata_subbasin(:,selsubbas,(sspidx(ssp):sspidx(ssp+1)-1))),c_ssp_main(ssp,:),3,monexceedances_xaxis);
            %scatter(1:12,squeeze(myFDC(:,selsubbas,(sspidx(ssp):sspidx(ssp+1)-1))),5,c_ssp_main(ssp,:),'o');
            ii=ii+1;
        end
        ylim([-56 206])
        applymyplotformat( strcat(charlbl(selsubbas),basinlabels.basinnames{selsubbas}))
        %add red hline at y=0
        yline(0,'Color','red')

        % Apply labels and formatting to graph
        xline(designQ_exceedances,'LineStyle',':','Color',mygraylines) % main design discharges
        %
        if selsubbas==1
            % add designQ labels
            topy=ylim;
            text(designQ_exceedances,floor(topy(2))*[1 1 1],compose("Q%d",designQ_exceedances),'HorizontalAlignment','center','FontWeight','bold')
            % add ylabels only in left plots
            ylabel("% change")
        end

        % add xlabels only in bottom plots
        xticks([0:25:100])
        if selsubbas<=4
            xticklabels([])
        end
        if selsubbas==6
            xlabel('Flow exceedance %')
        end

        if selsubbas~=1 && selsubbas~=5
            yticklabels([])
        end

        %decided that log scale didnt add more value
        %set(gca,'YScale','log')
    end
    legend(l3,rcpnames,'Orientation',"horizontal")

    sgtitle(sprintf("Historical and future FDC for %s",tf_full{seltf}))
end

%% Calculate min/max mean for each scenario for each subbasin
% For total pot
for selsubbas=1:8
    for selq=1:3
        idata=squeeze(designQ_selchanges_subbasin(selq,:,selsubbas))';
        rcpgroups=[0 repelem([1:6],4)]';

        % Save mean in 3d mat with rcp x pottype x subbasin
        idata_byrcp(:,:,selq)=[repelem(idata(1),4)' reshape(idata(2:end),[4  6])];
        idata_mean(:,selq,selsubbas)=groupsummary(idata,rcpgroups,'mean') ;
        idata_neg(:,selq,selsubbas)= groupsummary(idata,rcpgroups,'min') - idata_mean(:,selq,selsubbas);
        idata_pos(:,selq,selsubbas)= groupsummary(idata,rcpgroups,'max')- idata_mean(:,selq,selsubbas);
    end
    subbasdata_byrcp_prct{selsubbas}=idata_byrcp;

end
nrcpscens=size(idata_mean,1);
disp("Reshaped data for rcp vs corner")

%% FINAL: Range plot for % changes in sub-basin design discharges
subplottight = @(m,n,p) subtightplot (m, n, p, [0.15 0.02],[.15 .1],[.09 .02]); %[vgap hgap], marg_h -[low up],marg_w -[l r]

c_sel= cbrewer2("Set2",3); %cmap8_wong_gray;
xcenter=1:nrcpscens;
xshifted=[xcenter-.15;
    xcenter;
    xcenter+.15    ];
figure
for selsubbas=1:nsubbas
    subplottight(2,4,selsubbas)
    hold all

    for selq=1:3
        ee1(selq)=errorbar(xshifted(selq,:),idata_mean(:,selq,selsubbas), idata_neg(:,selq,selsubbas), idata_pos(:,selq,selsubbas),'_', ...
            'MarkerSize',0.1,'CapSize',0,'LineWidth',2,'Color',c_sel(selq,:),'DisplayName',designQnames{selq});
        %ee2=plot(xshifted(selq,:),idata_mean(:,selq,selsubbas),"_k",'DisplayName','Mean of corners'); %,'Color',mygraylines)
        ee2=plot(xshifted(selq,:),subbasdata_byrcp_prct{selsubbas}(:,:,selq),".",'Color',mygraylines,'DisplayName','4 corner models'); %,'Color',mygraylines)

    end
    xline([1.5 4.5]); %,'LineWidth',1.5)
    applymyplotformat( strcat(charlbl(selsubbas),basinlabels.basinnames{selsubbas}))
    set(gca,'XGrid','off','YMinorGrid','on' )
    %add red hline at y=0
    yline(0,'Color','red')

    if selsubbas<9
        ylim([-50 150])
        % ylim([0.1 1100]) % if theory pot incl
    end

    % add xlabels
    xticks([1:7])
    xlim([1.5 7.5])
    %if selsubbas<=4
    %   xticklabels([])
    %else
    xticklabels(["Historical" rcpnames rcpnames])
    % add mid-far label
    text([3 6] ,130*[1 1],tname,'HorizontalAlignment','center')%,'FontAngle','italic')
    %end
    % add ylabels only in left plots
    if selsubbas~=1 && selsubbas~=5
        yticklabels([])
    else
        ylabel("% change")
    end
end
legend([ee1 ee2(1)],'Location','southoutside', 'Orientation', 'horizontal')