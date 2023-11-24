% Try create spatial map w Onions for diff potential type

clear all
close all
addpath(genpath('G:\SurfDrive\GitConnect\Hydrus\devFiles\'), ...
    genpath('G:\SurfDrive\HPmodel\Hydrus\'))
run('myVarNames.m')
load('G:\SurfDrive\HPmodel\output\Figs_trial\MainScenarios.mat', 'basindata','catchments_cl',...
    'subPot_GWh', 'subPot_num', 'subPot_GWh_perarea','subPot_GWh_percapita')
load('G:\SurfDrive\HPmodel\data\ASIA\Basin_UIB\PantpeBasin_103.mat', 'outside')
catchments_cl=maskBasin(catchments_cl,~outside);
selsearch=3; % full mixed sust scenario
trialfigs=0;
cmap8=brighten(cmap8_wong,0);

label_subbasin_basin =[ basindata.basinnames{:},"ALL Indus"];

%% Manually assign xy for onions
%'Kabul'  'Swat'  'Indus''Jhelum''Chenab'    {'Ravi'}    {'Beas'}  {'Satluj'} 'All basins'
bas_x =[670     1050     1800    1360   1600    1400    1900    2300 2400];
bas_y =[410     380     500     600     890    1200    1100    1250 300] ;

%% FINAL: Onion plot for total potential in GWh/yr
for selsearch=1:3
    my_bubble_size=sqrt(subPot_GWh);
    selcols=1:4;
    figure;
    h=imagescnan(catchments_cl);hold all
    colormap(cmap8)
    set(h, 'AlphaData', subbasinalpha*(catchments_cl>0))
    %Add theoretical
    bubblechart(bas_x,bas_y,my_bubble_size(:,1,selsearch),cmap6_pottype(1,:),"MarkerFaceAlpha",0.6)
    %Add TFS
    for ii=selcols(2:end)
        hold on
        bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0.3)
    end
    ii=ii+1;
    %Add visualized
    bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0,'LineWidth',1.5)
    %Add existing
    ii=ii+1;
    bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0,'LineWidth',1.5,'DisplayName',"Existing")
    legend(pottypes6,'Location','northwest')
    xticklabels(basindata.basinnames)
    bubblesize([5, 200])
    blgd=bubblelegend('GWh/yr','Style','telescopic','Box','off','Location','southwest');
    title(sprintf("Sub-basin total potential for %s search",searchtypes{selsearch}))
    text(bas_x,bas_y,label_subbasin_basin,'HorizontalAlignment','center','FontAngle','italic')
    %colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'Ticklabels',(basindata.basinnames),'TickDirection','none')
    % Scale bubble size legend - need to do this manually!
    tmp=blgd.LimitLabels;
    blgdlims=[str2num(tmp{1}) str2num(tmp{2})];
    tmp2=string(compose("%.0f",blgdlims.^2));
    blgd.LimitLabels = tmp2;
end

%% FINAL: Onion plot for per capita potential in MWh/yr
for selsearch=1:3
    my_bubble_size=subPot_GWh_percapita*1000;
    selcols=1:4;
    figure;
    h=imagescnan(catchments_cl);hold all
    colormap(cmap8)
    set(h, 'AlphaData', subbasinalpha*(catchments_cl>0))
    %Add theoretical
    bubblechart(bas_x,bas_y,my_bubble_size(:,1,selsearch),cmap6_pottype(1,:),"MarkerFaceAlpha",0.6)
    %Add TFS
    for ii=selcols(2:end)
        hold on
        bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0.3)
    end
    ii=ii+1;
    %Add visualized
    bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0,'LineWidth',1.5)
    %Add existing
    ii=ii+1;
    bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0,'LineWidth',1.5,'DisplayName',"Existing")
    legend(pottypes6,'Location','northwest')
    xticklabels(basindata.basinnames)
    bubblesize([5, 200])
    blgd=bubblelegend("MWh/yr per capita",'Style','telescopic','Box','off','Location','southwest');
    title(sprintf("Sub-basin MWh/yr per capita for %s search",searchtypes{selsearch}))
    text(bas_x,bas_y,label_subbasin_basin,'HorizontalAlignment','center','FontAngle','italic')
    %colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'Ticklabels',(basindata.basinnames),'TickDirection','none')
end

%% FINAL: Onion plot for specific potential in GWh/yr per km2
for selsearch=3% 1:3
    my_bubble_size=subPot_GWh_perarea*1e6; %convert from per m2 to per km2
    selcols=1:4;
    figure;
    h=imagescnan(catchments_cl);hold all
    colormap(cmap8)
    set(h, 'AlphaData', subbasinalpha*(catchments_cl>0))
    %Add theoretical
    bubblechart(bas_x,bas_y,my_bubble_size(:,1,selsearch),cmap6_pottype(1,:),"MarkerFaceAlpha",0.6)
    %Add TFS
    for ii=selcols(2:end)
        hold on
        bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0.4)
    end
    ii=ii+1;
    %Add visualized
    bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0,'LineWidth',1.5)
    %Add existing
    ii=ii+1;
    bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0,'LineWidth',1.5,'DisplayName',"Existing")
    legend(pottypes6,'Location','northwest')
    xticklabels(basindata.basinnames)
    bubblesize([5, 200])
    blgd=bubblelegend("GWh/yr per km^2",'Style','telescopic','Box','off','Location','southwest');
    title(sprintf("Sub-basin GWh/yr per km^2 for %s search",searchtypes{selsearch}))
    text(bas_x,bas_y,label_subbasin_basin,'HorizontalAlignment','center','FontAngle','italic')
    %colorbar('Ticks', [101.4 102.2 103.1 104 105 105.8 106.7 107.6],'Ticklabels',(basindata.basinnames),'TickDirection','none')
end

%% Try onion w manual circle sizing
figure
selsearch=3;
my_bubble_size=subPot_GWh_percapita*1e3;
selcols=1:4;
subplot(2,2,1)
h=imagescnan(catchments_cl);hold all
colormap(cmap8)
set(h, 'AlphaData', subbasinalpha*(catchments_cl>0))
%Add theoretical
bubblechart(bas_x,bas_y,my_bubble_size(:,1,selsearch),cmap6_pottype(1,:),"MarkerFaceAlpha",0.6)
bubblesize
%Add TFS
for ii=selcols(2:end)
    hold on
    bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0.3)
    bubblesize
end
ii=ii+1;
%Add visualized
bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0)% ,'LineWidth',1.5)
bubblesize
%Add existing
ii=ii+1;
bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0)%,'LineWidth',1.5,'DisplayName',"Existing")
bubblesize
legend(pottypes6,'Location','northwest')
xticklabels(basindata.basinnames)
bubblesize([5, 120])
blgd=bubblelegend("MWh/yr per capita",'Style','telescopic','Box','off','Location','southwest');
sgtitle(sprintf("Sub-basin MWh/yr per capita for %s search",searchtypes{selsearch}))
title("Default size")
text(bas_x,bas_y,label_subbasin_basin,'HorizontalAlignment','center','FontAngle','italic')

subplot(2,2,2)
h=imagescnan(catchments_cl);hold all
colormap(cmap8)
set(h, 'AlphaData', subbasinalpha*(catchments_cl>0))
%Add theoretical
bubblechart(bas_x,bas_y,my_bubble_size(:,1,selsearch),cmap6_pottype(1,:),"MarkerFaceAlpha",0.6)
%Add TFS
for ii=selcols(2:end)
    hold on
    bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0.3)
end
ii=ii+1;
%Add visualized
bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0)% ,'LineWidth',1.5)
%Add existing
ii=ii+1;
bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0)%,'LineWidth',1.5,'DisplayName',"Existing")
%legend(pottypes6,'Location','northwest')
xticklabels(basindata.basinnames)
xx=[min(my_bubble_size(:,:,selsearch),[],'all') max(my_bubble_size(:,:,selsearch),[],'all')];
bubblesize(sqrt(xx)*10)
blgd=bubblelegend("MWh/yr per capita",'Style','telescopic','Box','off','Location','southwest');
text(bas_x,bas_y,label_subbasin_basin,'HorizontalAlignment','center','FontAngle','italic')
title("size based on area scaling")

%
subplot(2,2,3)
selsearch=3;
Xcenter=bas_x';
Ycenter=bas_y';
Csize=rescale(my_bubble_size(:,:,selsearch),30,600);
%Csize=my_bubble_size(:,:,selsearch);

Colmap=cbrewer2('Oranges',6); %lines;
Couter=Csize(:,1);

%    figure;
imagescnan(catchments_cl);hold all
colormap(cmap8_light)

arrayfun(@(x,y,diameter) rectangle('Position', [x-diameter/2, y-diameter/2, diameter, diameter], 'Curvature', [1 1], 'EdgeColor', 'w','FaceColor',Colmap(1,:)), Xcenter, Ycenter, Couter);
for ii=2:4 % Tech, Fin, Sust
    arrayfun(@(x,y,diameter,diameter_outer) rectangle('Position', [x-diameter/2, y+diameter_outer/2-diameter, diameter, diameter], 'Curvature', [1 1],'EdgeColor',  'w','FaceColor',Colmap(ii,:)), Xcenter, Ycenter, Csize(:,ii),Couter);
end
ii=ii+1;
%Add visualized
arrayfun(@(x,y,diameter,diameter_outer) rectangle('Position', [x-diameter/2, y+diameter_outer/2-diameter, diameter, diameter], 'Curvature', [1 1],'EdgeColor', Colmap(ii,:),'LineStyle','--','LineWidth',1.1), Xcenter, Ycenter, Csize(:,ii),Couter);
%Add existing
ii=ii+1;
arrayfun(@(x,y,diameter,diameter_outer) rectangle('Position', [x-diameter/2, y+diameter_outer/2-diameter, diameter, diameter], 'Curvature', [1 1],'EdgeColor', 'k','LineWidth',1.5), Xcenter, Ycenter, Csize(:,ii),Couter);
title(sprintf("Sub-basin MWh/yr per capital for %s search",searchtypes{selsearch}))
text(bas_x,bas_y,label_subbasin_basin,'HorizontalAlignment','center','FontAngle','italic')
title("symbol diameter based on value")

%
subplot(2,2,4)
Csize=sqrt(my_bubble_size(:,:,selsearch))*60;

Colmap=cbrewer2('Oranges',6); %lines;
Couter=Csize(:,1);

%    figure;
imagescnan(catchments_cl);hold all
colormap(cmap8_light)

arrayfun(@(x,y,diameter) rectangle('Position', [x-diameter/2, y-diameter/2, diameter, diameter], 'Curvature', [1 1], 'EdgeColor', 'w','FaceColor',Colmap(1,:)), Xcenter, Ycenter, Couter);
for ii=2:4 % Tech, Fin, Sust
    arrayfun(@(x,y,diameter,diameter_outer) rectangle('Position', [x-diameter/2, y+diameter_outer/2-diameter, diameter, diameter], 'Curvature', [1 1],'EdgeColor',  'w','FaceColor',Colmap(ii,:)), Xcenter, Ycenter, Csize(:,ii),Couter);
end
ii=ii+1;
%Add visualized
arrayfun(@(x,y,diameter,diameter_outer) rectangle('Position', [x-diameter/2, y+diameter_outer/2-diameter, diameter, diameter], 'Curvature', [1 1],'EdgeColor', Colmap(ii,:),'LineStyle','--','LineWidth',1.1), Xcenter, Ycenter, Csize(:,ii),Couter);
%Add existing
ii=ii+1;
arrayfun(@(x,y,diameter,diameter_outer) rectangle('Position', [x-diameter/2, y+diameter_outer/2-diameter, diameter, diameter], 'Curvature', [1 1],'EdgeColor', 'k'), Xcenter, Ycenter, Csize(:,ii),Couter);
title(sprintf("Sub-basin MWh/yr per capital for %s search",searchtypes{selsearch}))
text(bas_x,bas_y,label_subbasin_basin,'HorizontalAlignment','center','FontAngle','italic')
title("symbol size based on area scaling")



%% Try rectangle for total basin potential

%% Try different ways to create onion plot
if trialfigs
    %% Concentric bubbles on a straight line
    my_bubble_size=sqrt(subPot_GWh);
    selcols=1:4;
    cmap6_pottype=cbrewer2('Oranges',5); %lines;
    figure;
    bubblechart(1:8,ones(1,8),my_bubble_size(:,1,selsearch),cmap6_pottype(1,:))
    for ii=selcols(2:end)
        hold on
        bubblechart(1:8,ones(1,8),my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:))
    end
    legend(pottypes5(selcols))
    xticklabels(basindata.basinnames)
    bubblesize([10, 200])
    bubblelegend("GWh/yr",'Style','telescopic')
    title(sprintf("Sub-basin total potential for %s search",searchtypes{selsearch}))

    %% Concentric scatter on a map
    my_bubble_size=rescale(subPot_GWh,500,20000);
    %[bas_y, bas_x]= ind2sub(size(catchments_cl),basindata.outlet_idxs)	;

    selcols=1:4;
    cmap6_pottype=cbrewer2('Oranges',5); %
    cmap8= cmap8_light;
    figure;
    imagescnan(catchments_cl);hold all
    colormap(cmap8)
    h(1)=scatter(bas_x,bas_y,my_bubble_size(:,1,selsearch),cmap6_pottype(1,:),"filled","MarkerFaceAlpha",0.6); % plot theoretical
    for ii=selcols(2:end)
        hold on
        h(ii)= scatter(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"filled","MarkerFaceAlpha",0.6); %plot other classes
    end
    ii=ii+1;
    h(ii)= scatter(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0.6); %plot other classes

    legend(pottypes5)
    grid on
    title(sprintf("Sub-basin total potential for %s search",searchtypes{selsearch}))

    %xticklabels(basindata.basinnames)
    %bubblesize([10, 200])
    %bubblelegend("GWh/yr",'Style','telescopic')

    %% GOOD: Shifted scatter on a straight line
    my_bubble_size=rescale(subPot_GWh,1000,20000);
    my_bubble_radius=sqrt(my_bubble_size/pi);
    % eval bubble xys
    %my_y=[ones(8,1)*0 -my_bubble_radius(:,1,selpolscen) + my_bubble_radius(:,selcols(2:end),1)/2];
    my_y=[ones(8,1)*0 ones(8,1)*0 ones(8,1)*0 ones(8,1)*0];

    selcols=1:2;
    cmap6_pottype=cbrewer2('Oranges',5); %
    cmap6_pottype= lines;
    figure;
    %imagescnan(catchments_cl>0);hold all
    h(1)=scatter(1:8,my_y(:,1),my_bubble_size(:,1,selsearch),cmap6_pottype(1,:),"filled","MarkerFaceAlpha",0.6); % plot theoretical
    for ii=selcols(2:end)
        hold on
        h(ii)= scatter(1:8,my_y(:,ii),my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"filled","MarkerFaceAlpha",0.6); %plot other classes
    end
    legend(pottypes5(selcols))
    grid on
    %xticklabels(basindata.basinnames)
    %% GOOD: Rectangle based plots - cannot make the circles transparent
    %function plotOnion(xbottom,ybottom,csize,colmap)
    % XYbottom are column vectors.
    % Csize is a matrix w same num of rows as XY bottom. Cols of Csize need to
    % have data with successively declining vals
    selsearch=3;
    Xcenter=bas_x';
    Ycenter=bas_y';
    %Csize=(subPot_GWh_percapita(:,:,1)*1000);
    %Csize=sqrt(subPot_GWh_percapita(:,:,1)*1000);
    Csize=rescale(subPot_GWh_percapita(:,:,selsearch),30,600);

    %Csize=rescale(subPot_GWh(:,:,3),20,500);

    Colmap=cbrewer2('Oranges',6); %lines;
    Couter=Csize(:,1);

    figure;
    imagescnan(catchments_cl);hold all
    colormap(cmap8_light)

    arrayfun(@(x,y,diameter) rectangle('Position', [x-diameter/2, y-diameter/2, diameter, diameter], 'Curvature', [1 1], 'EdgeColor', 'w','FaceColor',Colmap(1,:)), Xcenter, Ycenter, Couter);
    for ii=2:4 % Tech, Fin, Sust
        arrayfun(@(x,y,diameter,diameter_outer) rectangle('Position', [x-diameter/2, y+diameter_outer/2-diameter, diameter, diameter], 'Curvature', [1 1],'EdgeColor',  'w','FaceColor',Colmap(ii,:)), Xcenter, Ycenter, Csize(:,ii),Couter);
    end
    ii=ii+1;
    %Add visualized
    arrayfun(@(x,y,diameter,diameter_outer) rectangle('Position', [x-diameter/2, y+diameter_outer/2-diameter, diameter, diameter], 'Curvature', [1 1],'EdgeColor', Colmap(ii,:),'LineStyle',':','LineWidth',1.5), Xcenter, Ycenter, Csize(:,ii),Couter);
    %Add existing
    ii=ii+1;
    arrayfun(@(x,y,diameter,diameter_outer) rectangle('Position', [x-diameter/2, y+diameter_outer/2-diameter, diameter, diameter], 'Curvature', [1 1],'EdgeColor', 'k','LineStyle','--','LineWidth',1.5), Xcenter, Ycenter, Csize(:,ii),Couter);
    title(sprintf("Sub-basin MWh/yr per capital for %s search",searchtypes{selsearch}))
    text(bas_x,bas_y,basindata.basinnames,'HorizontalAlignment','center','FontAngle','italic')

    %% Test data
    x = transpose(1:4);
    y = zeros(4,1);
    sz1=1000 + (1000-100).*rand(4,1),"filled";
    sz2=500 + (1000-50).*rand(4,1);
    sz3=50 + (500-5).*rand(4,1);

    figure;
    scatter(1:4, zeros(1,4),sz1 )
    hold all
    scatter(1:4, zeros(1,4),sz2 ,"filled")
    scatter(1:4, zeros(1,4), sz3,"filled")
    %
    figure
    viscircles([x y],sz1)

    %% Test w viscircles for GWh/yr per capital
    my_bubble_radius=rescale(subPot_GWh_percapita,10,300);
    fig=figure('Position', get(0, 'Screensize'));
    imagescnan(catchments_cl);hold all
    colormap(brighten(cmap8,.7))
    hold all
    for ii=selcols
        %hold on
        viscircles([bas_x' bas_y'],my_bubble_radius(:,ii,selsearch), 'Color',cmap6_pottype(ii,:), 'LineWidth',0.5); %,"MarkerFaceAlpha",0.3);
    end
    ii=ii+1;
    %Add visualized
    viscircles([bas_x' bas_y'],my_bubble_radius(:,ii,selsearch),'Color',cmap6_pottype(ii,:),'LineWidth',0.5,'LineStyle', ':');
    %Add existing
    ii=ii+1;
    viscircles([bas_x' bas_y'],my_bubble_radius(:,ii,selsearch),'Color','k','LineWidth',0.5);

    %% Try bubblesize
    my_bubble_size=subPot_GWh_percapita*1000;
    selcols=1:2;
    cmap6_pottype=cbrewer2('Oranges',5); %lines;
    figure;
    imagescnan(catchments_cl);hold all
    colormap(cmap8)
    bubblechart(bas_x,bas_y,my_bubble_size(:,1,selsearch),cmap6_pottype(1,:),"MarkerFaceAlpha",0.6)
    for ii=selcols(2:end)
        hold on
        bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap6_pottype(ii,:),"MarkerFaceAlpha",0.3)
    end
    % ii=ii+1;
    % bubblechart(bas_x,bas_y,my_bubble_size(:,ii,selsearch),cmap5_orange(ii,:),"MarkerFaceAlpha",0,'LineWidth',1.5)
    % legend(pottypes5)
    % xticklabels(basindata.basinnames)
    % bubblesize([10, 200])
    % blgd=bubblelegend("MWh/yr per capita",'Style','telescopic','Box','off','Location','southwest');
    % title(sprintf("Sub-basin MWh/yr per capital for %s search",searchtypes{selsearch}))
    % text(bas_x,bas_y,basindata.basinnames,'HorizontalAlignment','center','FontAngle','italic')

end

