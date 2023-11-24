% Fig_SubBasinPotential_Prctchange_allpots
% Created By    : Sanita Dhaubanjar
% Created For	: SustainIndus project WP2 Paper publication
%=========================
% Code for loading pre-processed output files from Hydrus model to create
% final plots used for publication
% Plot sub-basin wise total potential and % change in potential for far
% future

clc
clearvars
close all
addpath(genpath(fullfile(pwd,'Hydrus')))
run('myVarNames_Fut.m')

%% Load saved subbasin data for fut potentials
ofmatname=fullfile(rootof,"FutFigs_cl0","FutFigs_cl1","FigData",'MainFutScenarios.mat');
load(ofmatname,'subPot_GWh','potnames','catchments_cl')
disp("Loaded saved subPot data")
potnames_cl5=[potnames(1) newpots4cl];
%skip the sust-tech pot in excel
selpot=[1 2 3 5 6 7 ];

%% Calculate min/max mean for % change for each scenario for each subbasin
subPot_GWh_prct_change=(subPot_GWh-subPot_GWh(:,:,1))./subPot_GWh(:,:,1)*100;

% range for % change
for selsubbas=1:9
    pp=1; % because selpot is not serial
    for p=selpot(1:4)
        idata=squeeze(subPot_GWh_prct_change(selsubbas,p,:)); % basins are rowwise and scens are column wise - mid fut first, far fut second
        rcpgroups=[0 repelem([1:6],4)]';

        % Save mean in 3d mat with rcp x pottype x subbasin
        idata_byrcp_prct(:,:,pp)=[repelem(idata(1),4)' reshape(idata(2:end),[4  6])];
        idata_mean_prct(:,pp,selsubbas)=groupsummary(idata,rcpgroups,'mean') ;
        idata_neg_prct(:,pp,selsubbas)=groupsummary(idata,rcpgroups,'min') -idata_mean_prct(:,pp,selsubbas);
        idata_pos_prct(:,pp,selsubbas)=groupsummary(idata,rcpgroups,'max')- idata_mean_prct(:,pp,selsubbas);
        pp=pp+1;
    end
    subbasdata_byrcp_prct{selsubbas}=idata_byrcp_prct;
end

nscens=size(idata_mean_prct,1);
npots=size(idata_mean_prct(:,:,selsubbas),2);
disp("Reshaped data for rcp vs corner")

%% FINAL: Range plot for sub-basin % change in potentials - WITH theoretical
subplottight = @(m,n,p) subtightplot (m, n, p, [0.08 0.05],[.08 .08],[.08 .01]); %[vgap hgap], marg_h -[low up],marg_w -[l r]

c_sel= c_pot5; %hsv(4); %cmap8_wong_gray;
xcenter=1:7;
xshifted=[ xcenter-.30
    xcenter-.10
    xcenter+.10
    xcenter+.30];
figure
lblidx=1;
for selsubbas=[9 1:8]
    subplottight(3,3,lblidx)
    hold all

    for pp=1:npots
        ee1(pp)=errorbar(xshifted(pp,:),idata_mean_prct(:,pp,selsubbas), idata_neg_prct(:,pp,selsubbas), idata_pos_prct(:,pp,selsubbas),'_', ...
            'MarkerSize',0.1,'CapSize',0,'LineWidth',3.5,'Color',c_sel(pp,:),'DisplayName',pottypes6{pp});
        ee2=plot(xshifted(pp,:),idata_mean_prct(:,pp,selsubbas),"_k",'DisplayName','Mean of corners','MarkerSize',4 ); %,'Color',mygraylines)
    end
    set(gca,'xtick',1:7,'xticklabel',["Historical" rcpnames rcpnames]) %Xticklength
    xline([1.5 4.5],"Color",mygraylines); %,'LineWidth',1.5)
    applymyplotformat(strcat(charlbl(lblidx),label_subbasin_basin(selsubbas )))
    xlim([1.5 7.5]) % can skip the historical one
    set(gca,'XGrid','off','YMinorGrid','on' )

    %add red hline at y=0
    yline(0,'Color','red')

    %set diff ylim for all UIB and other subbasin plots
    if selsubbas<9
        ylim([-50 150])
        yticks(-50:25:150)
    end

    if ismember(lblidx,[1,4,7])
        ylabel("% change")
    end

    %     %keep yaxis only in left and remove from internal subplots
    %     if ~ismember(lblidx,[1,4,7])
    %         yticklabels([])
    %     else
    %         ylabel("% change")
    %     end
    % add annotations only to the UIB subplot
    if selsubbas==9
        set(gca,'YColor',extraaxis ,'XColor',extraaxis );
        % add mid-far label
        topy=ylim;%text(1,max(idataTWh(1,:)),strcat("\leftarrow ", tname(1),"   " ,tname(2)," \rightarrow "),'HorizontalAlignment','center')%,'FontAngle','italic')
        text(4.5,floor(topy(2)),strcat(tname(1),"      " ,tname(2)),'HorizontalAlignment','center')%,'FontAngle','italic')

    end
    %addAlphaLabel(lblidx,"inside")
    lblidx=lblidx+1;
end

legend([ee1 ee2(1)],'Location','southoutside','NumColumns', 3, 'Orientation', 'horizontal')
sgtitle("% change in potential")

%% Calculate min/max mean for TOTAL potential for each scenario for each subbasin - SKIP THEORETICAL
subPot_TWh=subPot_GWh/1000;

for selsubbas=1:9
    pp=1; % because selpot is not serial
    for p=selpot(1:4)
        idata=squeeze(subPot_TWh(selsubbas,p,:)); % basins are rowwise and scens are column wise - mid fut first, far fut second
        rcpgroups=[0 repelem([1:6],4)]';

        % Save mean in 3d mat with rcp x pottype x subbasin
        idata_byrcp_prct(:,:,pp)=[repelem(idata(1),4)' reshape(idata(2:end),[4  6])];
        idata_mean_prct(:,pp,selsubbas)=groupsummary(idata,rcpgroups,'mean') ;
        idata_neg_prct(:,pp,selsubbas)=groupsummary(idata,rcpgroups,'min') -idata_mean_prct(:,pp,selsubbas);
        idata_pos_prct(:,pp,selsubbas)=groupsummary(idata,rcpgroups,'max')- idata_mean_prct(:,pp,selsubbas);
        pp=pp+1;
    end
    subbasdata_byrcp_prct{selsubbas}=idata_byrcp_prct;
end

nscens=size(idata_mean_prct,1);
npots=size(idata_mean_prct(:,:,selsubbas),2);
disp("Reshaped data for rcp vs corner")
%% FINAL: Range plot for sub-basin TOTAL potentials - WITH theoretical
subplottight = @(m,n,p) subtightplot (m, n, p, [0.08 0.05],[.08 .08],[.08 .01]); %[vgap hgap], marg_h -[low up],marg_w -[l r]
skiptheorypot=2:4;
c_sel= c_pot5; %hsv(4); %cmap8_wong_gray;
xcenter=1:7;
xshifted=[ xcenter-.30
    xcenter-.10
    xcenter+.10
    xcenter+.30];
figure
lblidx=1;
for selsubbas=[9 1:8]
    subplottight(3,3,lblidx)
    hold all

    for pp=skiptheorypot
        ee1(pp)=errorbar(xshifted(pp,:),idata_mean_prct(:,pp,selsubbas), idata_neg_prct(:,pp,selsubbas), idata_pos_prct(:,pp,selsubbas),'_', ...
            'MarkerSize',0.1,'CapSize',0,'LineWidth',3.5,'Color',c_sel(pp,:),'DisplayName',pottypes6{pp});
        ee2=plot(xshifted(pp,:),idata_mean_prct(:,pp,selsubbas),"_k",'DisplayName','Ensemble mean','MarkerSize',4 ); %,'Color',mygraylines)
    end
    set(gca,'xtick',1:7,'xticklabel',["Historical" rcpnames rcpnames]) %Xticklength
    xline([1.5 4.5],"Color",mygraylines); %,'LineWidth',1.5)
    applymyplotformat(strcat(charlbl(lblidx),label_subbasin_basin(selsubbas )))
    xlim([1.5 7.5]) % can skip the historical one
    set(gca,'XGrid','off','YMinorGrid','on' )
  
    %set diff ylim for all UIB and other subbasin plots
    if selsubbas<9
        ylim([0 100])
    %    yticks(-50:25:150)
    end

    if ismember(lblidx,[1,4,7])
        ylabel(energylabel)
    end
    % add annotations only to the UIB subplot
    if selsubbas==9
        set(gca,'YColor',extraaxis ,'XColor',extraaxis );
        % add mid-far label
        topy=ylim;%text(1,max(idataTWh(1,:)),strcat("\leftarrow ", tname(1),"   " ,tname(2)," \rightarrow "),'HorizontalAlignment','center')%,'FontAngle','italic')
        text([3,6],floor(topy(2))*[1 1],tname,'HorizontalAlignment','center')%,'FontAngle','italic')

    end
    %addAlphaLabel(lblidx,"inside")
    lblidx=lblidx+1;
end

legend([ee1(skiptheorypot) ee2],'Location','southoutside', 'Orientation', 'horizontal')
%sgtitle("Total potential")
%%
disp("***************************************EOF***************************************")
