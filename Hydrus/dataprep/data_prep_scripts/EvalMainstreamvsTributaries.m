% Delineate final basin domain excluding endorheic part and define the
% mainstreams vs tributaries for use with hydrus
% Load the streamorder files prepped in GIS. Streamordering was tried w
% accChannel+Strahler, Qchannel+Strahler and Qchannel+Qstreamorder. The
% last option is used for final run

clear all
close all

% Load basin data
root= pwd;%"G:\SurfDrive\GitConnect";
matpath=fullfile(root,"data","UI","data");
dataout500m='UI500m_ArcGIS.mat'; % append to existing
matout=fullfile(matpath,dataout500m);
%
res="500m"; cellsz_m=500;
%
main_do_m= 4*1e3;
trib_do_m= 2*1e3; %segment length in m
%
ch_method=2;   % Channel identification method 1 = based on acc, 2=based on Q
getstrahler=1; % Arcgis based strahler for stream order =1, Qbased stream order=2  
trypointgen=0; % test outlet generation

load(fullfile(matpath,'UI500m_ArcGIS.mat'), 'acc','adir', 'outlets','outlet_orig','catchments','basinlabels', 'flowdist');
load(fullfile(matpath,'MLAB500m_40yrClimatology_m3day_R3.mat'), 'Q500m_m3day')
Qannual_m3s=Q500m_m3day{13}/(24*60*60);

%% Disable small basins
figure
subplot(1,2,1)
imagescnan(catchments)
%subbasinmask_uniont=subbasinmask_union;
catchments(catchments>=109)=0;
outside=catchments==0;
title("catchments")

subplot(1,2,2)
imagescnan(outside)
title("outside")
% Select only outlets that matter
idx=find(outlets & outlets <=109); %find only outlets that matter
ptID=outlets(idx);
[ir0,ic0]=ind2sub(size(outlets), idx);

%% Get channel cells
if ch_method==1
    [rch,cch,ich,channel]=getChannel(acc, outside);
    fprintf('Channel based on acc makes up %0.2f%% of basin cells\n',sum(channel(:))/sum(~outside(:))*100);
    fprintf('Q>=0.1m^3/s make up %0.2f%% of basin cells\n',sum((Qannual_m3s)>=0.1 & ~outside,'all')/sum(~outside(:))*100);
else
    channel =  Qannual_m3s>=0.1 & ~outside;
    ich=find(channel);
    [rch,cch]=ind2sub(size(acc), ich);
    fprintf('Channel based on Q>=0.1m^3/2 makes up %0.2f%% of basin cells\n',sum(channel(:))/sum(~outside(:))*100);
end

%% Find mainstream for all selected subbasin outlets
% get r,c for points on the mainstream and for
for pt=1:length(ptID)
    inlet=idx(pt);
    mainstream = find_mainstream(Qannual_m3s,inlet,adir);
    [r0,c0]=ind2sub(size(acc), mainstream);
    
    % Store R,C
    if pt==1
        main_c= c0;
        main_r= r0;
        main_i= mainstream;
    else
        main_c=[main_c; c0];
        main_r=[main_r; r0];
        main_i=[main_i; mainstream];
    end
    
end

%% need to potentially add a step where we check for overlap betn channel
% and mainstream cells
% tmp=zeros(size(channel));
% tmp(main_i)=1;
% figure;imagescnan(tmp)

% Only take main_i that lies in the channel, i.e., has Q >=0.1m3/s
main_i_valid= main_i(channel(main_i)==1 & outside(main_i)==0);
%% Define main and trib channels
channel_main_trib = single(channel); % tributaries=1
channel_main_trib(main_i_valid)=2;  % mainstream=2

countUniques(channel_main_trib(:))
if all(channel_main_trib>0==channel)
    disp("Channel and main vs trib match")
else
    disp("!!!!Channel and main vs trib DONT match!!!!")
end
    % figure;imagescnan(channel_main_trib)
% title("Mainstream and tributary")

%% FINAL: plot subcatchments with mainstream and tributaries and subbasin outlets
cmap=cbrewer('qual', 'Set3', length(ptID)-1);
catchments(catchments==0)=nan;
%
figure %('Position', get(0, 'Screensize'))
imagescnan(catchments)
hold all;
% plot channels, mainstream and outlets
%l1= plot(ic0,ir0,"ok",'markersize',15) ;  %plot outlets
l2= plot(cch,rch,'.b','markersize',3);%,'Color', '#89D1FE') ;  %plot channel in light blue
l3= plot(main_c,main_r,'.r','markersize',8);%,'Color', [0.5 0 0]) ; %plot mainstream - skip the 0
colormap(cmap)
% Add basin names to clabel
c=colorbar('Ticks', [basinlabels.basinIDs(1:8)],'TickLabels',{basinlabels.basinnames{1:8}},'Direction','reverse');
c.Label.String="Sub-basins";
c.TickLength=0;
xticks([]); yticks([]);
legend([l3 l2],'Mainstreams','Tributaries')
%export_fig('G:\SurfDrive\GitConnect\data\UI\MainstreamVsTributary.pdf')

%% Generate and plot outlets
if trypointgen
    [r1,c1] = find_outletsInChannel(channel_main_trib==1,flowdist,trib_do_m/cellsz_m);
    [r2,c2] = find_outletsInChannel(channel_main_trib==2,flowdist,main_do_m/cellsz_m);
    
    % add outlets to previous plots;
    hold on;
    l4=plot(c2,r2, 'xk','markersize',4);
    l5=plot(c1,r1, '.r','markersize',4);% ,'Color',.1*[1 1 1]);
    l6=plot(0,0, '.r','markersize',8);% ,'Color',.1*[1 1 1]);
    
    legend([l1 l3 l4 l2 l5],'Outlets', 'Mainstream Channel', 'Mainstream Outlet', 'Tributary Channel', 'Tributary Outlet')
    %    legend([l3 l4 l2 l6], 'Mainstream', 'Mainstream Outlet', 'Tributary', 'Tributary Outlet')
end

%% Process streamorder
if getstrahler
    if ch_method==1 %acc based channel and strahler
        load(fullfile(matpath,'streamorder_STRAHLER.mat'), 'data');
        streamorder_STRAHLER=data;
    elseif getstrahler==1   %Q based channel and strahler
        load(fullfile(matpath,'streamorder_STRAHLER_Qbased.mat'), 'data');
        streamorder_STRAHLER=data;
    % map order = 1 to tertiary(1), order = 2 to secondary(2) and order =3,4,5% to primary(3)
    channel_ord=changem(streamorder_STRAHLER,[1 1 2 2 3 3 3], 1:7 );
    countUniques(channel_ord)
    elseif getstrahler==2 % Qbased streamorder
        load(fullfile(matpath,'streamorder_Qclasses.mat'), 'data');
        streamorder_STRAHLER=data;
    end

    figure;
    imagescnan(channel_ord)
    fprintf('Channel makes up %0.2f%% of basin cells\n',sum(channel_ord(:)>0)/sum(~outside(:))*100);
    
    %% FINAL: plot subcatchments with streamorder and subbasin outlets
    [r_ter, c_ter]=find(channel_ord==1);
    [r_sec, c_sec]=find(channel_ord==2);
    [r_pri, c_pri]=find(channel_ord==3);
    
    m_sz=5;
    cmap=cbrewer('qual', 'Pastel2', length(ptID)-1);
    catchments(catchments==0)=nan;
    %
    figure%('Position', get(0, 'Screensize'))
    imagescnan(catchments)
    hold all;
    % plot channels, mainstream and outlets
    l1= plot(c_pri,r_pri,".",'markersize',m_sz,'DisplayName','Primary') ;
    l2= plot(c_sec,r_sec,'.','markersize',m_sz,'DisplayName','Secondary');%,'Color', '#89D1FE') ;  %plot channel in light blue
    l3= plot(c_ter,r_ter,'.','markersize',m_sz-2,'DisplayName','Tertiary');%,'Color', [0.5 0 0]) ; %plot mainstream - skip the 0
    colormap(cmap)
    % Add basin names to clabel
    c=colorbar('Ticks', [basinlabels.basinIDs(1:8)],'TickLabels',{basinlabels.basinnames{1:8}},'Direction','reverse');
    c.Label.String="Sub-basins";
    c.TickLength=0;
    xticks([]); yticks([]);
    legend()
        %set(gca,'FontSize',14,'Layer','top');
    export_fig(fullfile(matpath,'\StreamOrder_QclassesNotStrahler.pdf'))
end
%% Save to tif
%savemat2Pantpetiff('G:\SurfDrive\GitConnect\data\data_prep\500m\channel_main_trib_500m_UIprj.tif',channel_main_trib)
% savemat2Pantpetiff('G:\SurfDrive\GitConnect\data\data_prep\500m\channel_main_trib_500m_Qbased_UIprj.tif',channel_main_trib)
% savemat2Pantpetiff('G:\SurfDrive\GitConnect\data\data_prep\500m\channel_Qbased_UIprj.tif',channel)

%% Save mainstream to existing UIB file
if ch_method==1
    save(matout,'channel','outside','channel_main_trib','channel_ord', '-append')
else
    matout=fullfile(matpath,"channel_Qbased.mat");
    save(matout,'channel','channel_main_trib','channel_ord')
end
disp('Processed data is saved')
disp("***************************************EOF***************************************")