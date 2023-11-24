function plotCostCurveScatter(COEAlls, PnetAlls, SysIDs,mycolor,myfill,mylinetype,linelabel)
%% Create a cost curve distinguishing between the two project types
% and scaling symbol size by PnetAlls

Techlim=0.5;% $/kWh
Econlim=0.1;% $/kWh

% Sort all cells in the order of the COEs from low to high
[COEISort,ind]=sort(COEAlls);
SysIDsSort=SysIDs(ind);
PnetSort=PnetAlls(ind); %in GWh
PnetISort_cum=cumsum(PnetSort);
% PnetISort_cum=[0; PnetISort_cum];
% COEISort=[0; COEISort];

%create symbol size vector
PnetSortSize= rescale(PnetSort,25,1500); %,'InputMax',maxscal,'InputMin',minscal);

%create color
if exist('mycolor')
    mycolor1=mycolor;
    mycolor2=mycolor;
    mycolor3=mycolor;
else
    mycolor1 =[0.8 0 0];%myred
    mycolor2=[0 0.5 0];%mygreen
    mycolor3=[0 0.45 0.74];    %riverblue
    figure;
end

if ~exist('mylinetype')
    mylinetype='-';
end
% GOOD: Cost curve simple
hold all
% add points for DP and RP
%plot(PnetISort_cum(SysIDsSort==2)/1000,COEISort(SysIDsSort==2),'o','Color',mycolor2,'DisplayName','River power');%'--','LineWidth',2,
%plot(PnetISort_cum(SysIDsSort==1)/1000,COEISort(SysIDsSort==1),'x','Color',mycolor1,'DisplayName','Diversion');%'--','LineWidth',2,
if myfill
    scatter(PnetISort_cum(SysIDsSort==2)/1000,COEISort(SysIDsSort==2),PnetSortSize(SysIDsSort==2),'^','filled','MarkerFaceAlpha',0.3,'MarkerFaceColor',mycolor2,'DisplayName','River power');%'--','LineWidth',2,
    scatter(PnetISort_cum(SysIDsSort==1)/1000,COEISort(SysIDsSort==1),PnetSortSize(SysIDsSort==1),'o','filled','MarkerFaceAlpha',0.3,'MarkerFaceColor',mycolor1,'DisplayName','Diversion');%'--','LineWidth',2,
else
    scatter(PnetISort_cum(SysIDsSort==2)/1000,COEISort(SysIDsSort==2),PnetSortSize(SysIDsSort==2),'^','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',mycolor2,'DisplayName','River power','LineWidth',1.15);%'--','LineWidth',2,
    scatter(PnetISort_cum(SysIDsSort==1)/1000,COEISort(SysIDsSort==1),PnetSortSize(SysIDsSort==1),'o','MarkerEdgeAlpha',0.3,'MarkerEdgeColor',mycolor1,'DisplayName','Diversion','LineWidth',1.15);%'--','LineWidth',2,
end
stairs(PnetISort_cum/1000,COEISort,mylinetype,'LineWidth',1.75,'Color',mycolor3,'DisplayName','Cost curve' );

if exist('linelabel','var')
    text(max(PnetISort_cum)/1000,max(COEISort),sprintf(' %s',linelabel),'Color',mycolor3,'FontAngle', 'italic','Rotation',90)%,'VerticalAlignment','baseline')
end

if ~exist('mycolor')
    % Add lines for potential cutoffs
    yline(Techlim,'LineStyle','--','DisplayName','Technical potential limit'); %Tech potential
    yline(Econlim,'LineStyle','--','Color','r','DisplayName','Economic potential limit'); %Economic potential
    %text(10,Econlim*1.1, "\it Economic potential limit")

    TotPotential = sum(PnetAlls(~isnan(PnetAlls)));
    TechPotential = sum(PnetAlls(COEAlls<Techlim));
    EconPotential =sum(PnetAlls(COEAlls<Econlim));
    maxpot = sprintf('Total potential: %0.1f TWh \nTechnical potential: %0.1f TWh \nEconomic potential: %0.1f TWh',TotPotential/1000,TechPotential/1000,EconPotential/1000);
    title(maxpot, 'FontWeight','bold','FontSize',10)
    ylabel('Per unit production cost (USD2010 per KWh)','FontSize',10)
    xlabel('Cumulative Energy (TWh)','FontSize',10)
    set(gca,'FontSize',10,'YScale','log')
    grid on; box on
    legend('-DynamicLegend','location','Best')
end
end
