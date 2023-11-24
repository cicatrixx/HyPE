function lh=plotCostCurve(COEAlls, PnetAlls, SysIDs,mycolor,mylinetype,myalpha,linelabel,nopoints)
%% Create a cost curve distinguishing between the two project types.
%% divides PnetAlls/1000 so unit changed in plots
%RP=triangle DP=circle
% no symbol size scaling
Techlim=0.5;% $/kWh
Econlim=0.1;% $/kWh
[COEISort,ind]=sort(COEAlls);
SysIDsSort=SysIDs(ind);

PnetISort=PnetAlls(ind);
TotPotential = sum(PnetAlls(~isnan(PnetAlls)));
TechPotential = sum(PnetAlls(COEAlls<Techlim));

EconPotential =sum(PnetAlls(COEAlls<Econlim));
PnetISort_cum=cumsum(PnetISort,'omitnan');
% PnetISort_cum=[0; PnetISort_cum];
% COEISort=[0; COEISort];

%create color
if exist('mycolor')
    mycolor1=mycolor;
    mycolor2=mycolor;
else
    mycolor1 =[0.8 0 0];%myred
    mycolor2=[0 0.5 0];%mygreen
    %riverblue=[0 0.45 0.74];
    figure;
end

if ~exist('mylinetype')
    mylinetype='-';
end

if ~exist('myalpha')
    myalpha=0.4;
end
% GOOD: Cost curve simple
hold all
% add points for DP and RP

%plot(PnetISort_cum(SysIDsSort==2)/1000,COEISort(SysIDsSort==2),'^','Color',mycolor2,'DisplayName','River power');%'--','LineWidth',2,
%plot(PnetISort_cum(SysIDsSort==1)/1000,COEISort(SysIDsSort==1),'o','Color',mycolor1,'DisplayName','Diversion');%'--','LineWidth',2,
if ~nopoints
scatter(PnetISort_cum(SysIDsSort==2)/1000,COEISort(SysIDsSort==2),'^','MarkerEdgeColor',mycolor2,'DisplayName','River power');%'--','LineWidth',2,
scatter(PnetISort_cum(SysIDsSort==1)/1000,COEISort(SysIDsSort==1),'o','MarkerEdgeColor',mycolor1,'DisplayName','Diversion');%'--','LineWidth',2,
end
lh=stairs(PnetISort_cum/1000,COEISort,mylinetype,'LineWidth',1.25,'Color',mycolor,'DisplayName','Cost curve' );
alpha(myalpha)

if exist('linelabel','var')
    text(max(PnetISort_cum)/1000-5,max(COEISort)-8,sprintf(' %s',linelabel),'Color',mycolor1,'FontAngle', 'italic','Rotation',90,'Interpreter','none'); %,'HorizontalAlignment','center')%,'VerticalAlignment','baseline')
end

if ~exist('mycolor')
    % Add lines for potential cutoffs
    yline(Techlim,'LineStyle','--','DisplayName','Technical potential limit'); %Tech potential
    yline(Econlim,'LineStyle','--','Color','r','DisplayName','Economic potential limit'); %Economic potential
    %text(10,Econlim*1.1, "\it Economic potential limit")

    maxpot = sprintf('Total potential: %0.1f TWh \nTechnical potential: %0.1f TWh \nEconomic potential: %0.1f TWh',TotPotential/1000,TechPotential/1000,EconPotential/1000);
    title(maxpot, 'FontWeight','bold','FontSize',10)
    ylabel('Per unit production cost (USD2010 per KWh)','FontSize',10)
    xlabel('Cumulative Energy (TWh)','FontSize',10)
    set(gca,'FontSize',10,'YScale','log')
    grid on; box on
    legend('-DynamicLegend','location','Best')
end
end
