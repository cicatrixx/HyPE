function [COETotminR, OptPTnet_GWh, OptPnet_kWh, OptP_MW, OptDH, OptDL, CostElements, OptInv, Pmin_W, OptPop, OptLV, OptSpecCap] = RDcostmodel(dfldamRD,dfhdamRD,dfPopCost,Qoutlets_design,Qoutlets_design_LF,Dis2Trans,LandValue,RDDepth,nbasin,outlet,cost_constr,cost_lim,Hazardrate_RD_prct,Dis2Road,makechanges2cost)
% Eval cost for all user provided combinations of DH and DL and identify
% the cheapest option.

% all outputs are for cost optimal dam sizing 
% COETotminR= COE
% OptPnet= Annual production in kWh of cost optimal dam 
% OptP= Plant capacity in MW
% OptDH= Dam height in m
% OptDL= Dam length in m
% CostElements= 14 COE components
% OptInv= Inv cost in $
% Pmin= Plant capacity in W
% OptPop= PopCost in $
% OptLV= LandValue cost in $
% OptSpecCap= specific investment cost $/kW
%
% b=k;
% dfldamRD=dfldamRD{b};
% dfhdamRD=single(dfhdamRD{b});
% dfPopCost=PopCost{b};
% Qoutlets_design=Qoutlets_design(b);
% Qoutlets_design_LF=Qoutlets_design_LF(b);
% Dis=DisOutlet(b);
% LandValue=LandValueRDlake{b};
% RDDepth=RDDepth(b);
% nbasin;
% outlet=b;
% cost_constr;
% cost_lim;
% Quakerate=0;

% Veileder 2012 Small hydropower http://publikasjoner.nve.no/veileder/2012/veileder2012_02.pdf

%% Variables from outer space
Qtot    = Qoutlets_design;      %Q_design based on monthly Q
LF      = Qoutlets_design_LF;   %Loadfactor based on monthly pattern and Qdesign
% dfldamRD,dfhdamRD are vectors
DH      = dfhdamRD;             %Dam height (m)
DL      = dfldamRD;             %Dam length (m)
PopCost = dfPopCost;            %Lake cost (GDP * ResettlementPop) + 10% for Benefit Sharin
%RD_minH=4; % minimum head in m
Hazardrate_RD=Hazardrate_RD_prct/100;
Hgross = DH-RDDepth; %Gross head vector: Dam height minus riverdepth to determine actual pressurized head

%% RDcostmodel constants
% for i=1:numel(DH) %not needed a minH check already done in main code
%     if Hgross(i)<RD_minH; Hgross(i)=0; end; % Minimum dam height of 3m
% end
run('config_costmodels.m')
% Make changes to config vals if any
if exist('makechanges2cost','var')
    %disp("Changes made to costconfig")
    eval(makechanges2cost);
end
%% Physical calculations
QQ=Qtot;
% Vectors w all possible options
%PTgross_GWh = rho*g*Qdesign_mean*Hgross*8760*1e-9;    %Theoretical gross potential
PTnet_GWh = rho*g*Qtot*Hgross*8760*1e-9;      % Net potential
P_W = eta_generation_RD*eta_transmission_RD*rho*g*Hgross.*QQ;                 % Capacity based on calculated Q (W)
Pnet_kWh = (P_W*1e-3*8760).*LF;              % Yearly production with LF (kWh)

%% Turbine based on IRENA hydropower, empirical
%Different cost metrics based on capacity, Q, head and dam height
%CostGunduz = -0.0254 + 1.81*QQ + 0.0848*(Hgross-Ahf) + 0.0553*L;    %Cost based on Gunduz & Sahin No transmission and 100y flood Q
%CostORNL = 1.2 * 110168 * DH^-0.35*(P*1e-3).^-0.3;                  %Initial capital cost ORNL formula ($/kW)
TurIrena = 1e6 *(1.1943* (P_W*1e-6).^0.7634) * IR2005;                          %Turbine investment costs ($) IRENA hydropower, empirical R2=0.94
% Pelton2jV = 1105.2*QQ.^-0.5106;                                       %Pelton turbine Veileder p156 2-jetter high head Q<10 h=1000
% Pelton6jV = (1559.6*QQ.^-0.5179).*(P*1e-3) ;                          %Pelton turbine Veileder p156 6-jetter high head 10>Q<35 h=1000 ($/kW)
% FrancisV = (1439.1*QQ.^-0.3143).*(P*1e-3);                                        %Francis turbine Veileder p157 0>Q<160 h=300 ($/kW)
% KaplanV = 11730.5*QQ.^-0.2953;                                        %Kaplan turbine Veileder p158 0>Q<400 h=15 ($/kW)
% SpecTurcost = TurIrena./(P*1e-3);                                    %Specific investment costs ($/kW)

%% Dam cost (RCC) based on Veileder report p62
%DV = DL * 0.5*(DW*DH); %Dam volume (m3) 6e6 m3 (assumed it is a triangle)
DCn = 0.72*DH.^1.8;     %Price (1000 NOK/m) RCC >1M m3
%DCn = 1.69*DH.^1.68;   %Price (1000 NOK/m) RCC >0.1M m3
TDC = DCn.*DL*1e3*ER_noktoUSD2010;   %Price ($)
TDC(TDC<20e6) =20e6;    %Minimum startup cost
TDamkW = TDC./(P_W*1e-3); %Dam cost per kW ($/kW)

%% Underground powerstation (depends on blasted volume) Veileder p101
BlastedV_m2= 78*Hgross.^0.5.*QQ.^0.7.*npowerunits.^0.1;
BlastCost=BlastedV_m2*2250*ER_noktoUSD2010;  % $

%Overground powerstation cost
PSSnok = -0.0006*QQ.^2 + 0.67.*QQ - 6.95;    %mil NOK/m (head 10-40m, 1 powerunit)
PSSnok(PSSnok<20) = 20;                     % Not lower than 20mil NOK
PSS = PSSnok*ER_noktoUSD2010*1e6+BlastCost;                        % $
PSSkW = PSS./(P_W*1e-3);                      % $/kW

%% Total electrical mechanical costs based on Veileder p147 (2 power units)
PEM = (7.5748*(P_W*1e-6).^0.65618) * ER_noktoUSD2010 *1e6;       %Total electro-technical equipment $
PEMkW = PEM./(P_W*1e-3);                           %Total electro-technical equipment $/kW

%% Penstock cost Veileder p94 (surface penstock and steel pipes)
npipe_RD=(DL/65);   %Three Gorges has a pipe every 65m.
PENnok = (6*D_RD+9.5)*DH.*npipe_RD; %1000 NOK/m Veileder Eq. x Dam height x Number of pipes
PEN = PENnok*ER_noktoUSD2010*1e3;            %$
PENkW = PEN./(P_W*1e-3);          %$/kW

%% Fish pasage mitigation cost DOEwater p32: title: Estimation of Economic Parameters of U.S. Hydropower Resources
FishCost = 1.3e6 *(P_W*1e-6).^0.56 * IR2002;

%% Miscellaneous
Misc = (-38.795.*log(QQ) + 309.89).*(P_W*1e-3)*ER_noktoUSD2010; % $

%% Add provincial water use to misc costs
Misc = Misc + PWC_pak * ER_paktoUSD * Pnet_kWh;

%% DisCost - Connection to nearest transmission line + per unit transmission cost
DisCostNOK = Powerline_allocator(Dis2Trans,P_W*1e-6); % NOK
DisCost = DisCostNOK*ER_noktoUSD2010 + COETransCost*Pnet_kWh; % $2010

%% RoadCost - Connection to nearest road
RoadCostNOK = road_allocator(Dis2Road,P_W*1e-6); % NOK 
RoadCost =RoadCostNOK*ER_noktoUSD2010; % $2010

%% Calculate annulized costs
AnnFac = interest/(1-((1+interest)^(-lifetime_RD)));   % Annuity factor using interest (10%) and economic lifetime (40 years) = 10,23%
AnnualCost = AnnFac * (TurIrena + PEM + PSS + TDC + PEN + FishCost + Misc + PopCost + DisCost + RoadCost + LandValue);
%COE = AnnualCost./Pnet_kWh;                           %Annual Cost of electricity ($/kWh)

%Collecting cost information
AnTurR = AnnFac * TurIrena;     %Turbine (annualized $)
AnDamR = AnnFac * TDC;          %Dam (annualized $)
AnPenR = AnnFac * PEN;          %Penstock (annualized $)
AnElecR = AnnFac * PEM;         %Electro-technical equipment ($)
AnPSR = AnnFac * PSS;           %Powerstation($)
AnFish = AnnFac * FishCost;     %Fish Cost ($)
AnMisc = AnnFac * Misc;         %misc cost ($)
AnPopR = AnnFac * PopCost;      %Pop Cost ($)
AnRoad = AnnFac * RoadCost;     %Pop Cost ($)
AnDis = AnnFac * DisCost;       %Distance2dem cost
AnLandVal = AnnFac * LandValue;  %Land value Cost ($)
AnOMR = AnnFac* ((TurIrena+TDC+PEN+PEM+PSS+FishCost+Misc+PopCost+DisCost+RoadCost +LandValue) * OMshare_RD);         %OM cost as fraction of total investment costs
AnCtotR = AnTurR + AnDamR + AnPenR + AnElecR + AnPSR + AnFish + AnMisc + AnPopR + AnDis + AnLandVal + AnOMR;
AnCOwner = AnCtotR * Ownersrate_RD;
AnCHazard = AnCtotR * Hazardrate_RD;
AnCtotR = AnCtotR + AnCOwner + AnCHazard;

%% Calculate COE cost components
COETurR = AnTurR./Pnet_kWh;
COEDamR = AnDamR./Pnet_kWh;
COEPenR = AnPenR./Pnet_kWh;
COEElecR = AnElecR./Pnet_kWh;
COEPSR = AnPSR./Pnet_kWh;
COEFish = AnFish./Pnet_kWh;
COEMisc = AnMisc./Pnet_kWh;
COEPop = AnPopR./Pnet_kWh;
COEDis = AnDis./Pnet_kWh;
COERoad = AnRoad./Pnet_kWh;
COELandVal = AnLandVal./Pnet_kWh;

COEOM = AnOMR./Pnet_kWh;
COEOwner = AnCOwner./Pnet_kWh;
COEHazard = AnCHazard./Pnet_kWh;
COETotR = COETurR+COEDamR+COEPenR+COEElecR+COEPSR+COEFish+COEMisc+COEPop+COEDis+COERoad+COELandVal+COEOM+COEOwner+COEHazard; %COETransCost;

%% Nan 0s for min search for optimal dam height
for jj=1:numel(COETotR)
    if COETotR(jj)==0; COETotR(jj)=NaN;end
end

%% Find cost optimal dam
[c,b] =min(COETotR);
COETotminR = c;
OptDH = DH(b);
OptDL = DL(b);
OptP_MW = P_W(b)*1e-6; %MW
OptPnet_kWh = Pnet_kWh(b); % kWh
OptPTnet_GWh= PTnet_GWh(b); %GWh
OptInvAnn = AnCtotR(b); % $ Annualized investments
Pmin_W = P_W(b); % W
OptPop = PopCost(b); % $
OptLV = LandValue(b); % $

OptInv1 = TurIrena(b) + TDC(b) + PEN(b) + PEM(b) + PSS(b) ...
    + FishCost(b) + Misc(b) + PopCost(b) + DisCost(b) +RoadCost(b) + LandValue(b); % $ absolute investments
OptInvOM = OMshare_RD * OptInv1;
OptInvsOwners = (OptInv1 + OptInvOM) * Ownersrate_RD;
OptInvsHazard = (OptInv1 + OptInvOM) * Hazardrate_RD;
OptInv = OptInv1 + OptInvOM + OptInvsOwners + OptInvsHazard;

OptSpecCap = OptInv/(Pmin_W*1e-3); % $/kW

CostElements{1}=COETurR;
CostElements{2}=COEDamR;
CostElements{3}=COEOM;
CostElements{4}=COEPenR;
CostElements{5}=COEElecR;
CostElements{6}=COEPSR;
CostElements{7}=COEFish;
CostElements{8}=COEMisc;
CostElements{9}=COEPop;
CostElements{10}=COEDis;
CostElements{11}=COERoad;
CostElements{12}=COELandVal;
CostElements{13}=COEOwner;
CostElements{14}=COEHazard;

%arcCostElements=cat(1,CostElements{:})';

if cost_constr==1
    if COETotminR>cost_lim
        COETotminR=NaN; %Select river dam lower than x$/kWh
        OptPnet_kWh=NaN;
    else
    end
%     
%     h=figure(1);clf;
%     bar(1:numel(DH),[COETurR' COEDamR' COEOMR' COEPenR' COEElecR' COEPSR' COEFish' COEMisc' COEPop' COEDis' COELandVal' COEOwner'],0.5,'stack');
%     legend('Turbine', 'Dam', 'OM','Penstock','Electro','Powerstation','Fish','Misc','Pop','Distance','Land','Owners')
%     title('COE stacked');
%     ylabel('$/kWh');
%     xlabel('Dam height (m)');
%     axis([0,300,0,0.3])
%     
%     fprintf('Cost_constrain=on %d Outlet #%d, Cost = %.2f $/kWh, Cap = %.2f MW, DH = %.0f m\n',nbasin,outlet,COETotminR,OptP,OptDH)
%     
else
%     fprintf('Cost_constrain=off %d Outlet #%d, Cost = %.2f $/kWh, Cap = %.2f MW, DH = %.0f m\n',nbasin,outlet,COETotminR,OptP,OptDH)
    
end

%%
% if COETotminR>2
%     COETotminR=NaN; %Select river dam lower than 0.5$/kWh
%     OptPnet=NaN;

% elseif COETotminR<0.2
%     h=figure('Visible','off');
%     bar(1:numel(DH),[COETurR' COEDamR' COEOMR' COEPenR' COEElecR' COEPSR' COEFish' COEMisc' COELandVal' COEDis'],0.5,'stack');
%     legend('Turbine', 'Dam', 'OM','Penstock','Electro','Powerstation','Fish','Misc','Lake','Distance')
%     title('COE stacked');
%     ylabel('$/kWh');
%     xlabel('Dam height (m)');
%     axis([0,inf,0,1])
%
%     pathname = fileparts('Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\NAM\output\RD\');
%     figfile = fullfile(pathname, sprintf('RDam_costperkWh_b%d_outlet%d.png',nbasin,outlet));
%     saveas(h,figfile);

%     hh=figure('Visible','off');
%     bar(1:numel(DH),[AnTurR' AnDamR' AnOMR' AnPenR' AnElecR' AnLandVal'],0.5,'stack');
%     legend('Turbine', 'Dam', 'OM','Penstock','Electro','Lake')
%     title('Investments stacked');
%     ylabel('Annualized investments');
%     xlabel('Dam width *90m');
%
%     pathname = fileparts('Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\NAM\output\RD\');
%     figfile = fullfile(pathname, sprintf('RDam_Invest_b%d_outlet%d.png',nbasin,outlet));
%     saveas(hh,figfile);

%     fprintf('Cost of outlet #%d = %.2f $/kWh\n',outlet, COETotminR)
%     fprintf('Capacity of outlet #%d = %.2f MW\n',outlet, OptP)
%     fprintf('DH of outlet #%d = %.0f m\n',outlet, OptDH)
%     fprintf('DL of outlet #%d = %.0f m\n',outlet, OptDL)
%     fprintf('Q of outlet #%d = %.2f m3/s \n',outlet, Qtot)
%     fprintf('LF of outlet #%d = %.2f \n\n',outlet, LF)
%     fprintf('Pnet of outlet #%d = %.0f MWh\n',outlet, OptPnet*1e-3)
% else
%
% end;

% COETotminRIdx = find(isnan(COETotminR)); %Find NaN index
% OptPnet(COETotminRIdx)=NaN;           %Remove Pnet with COE=NaN


%                 clf(figure(1),'reset');
%                 clf(figure(2),'reset');
%                 clf(figure(3),'reset');
%                 clf(figure(4),'reset');
%
%                 figure(1);
%                 plot(DH,COETot);
%                 title('COETot');
%                 xlabel('Dam width (m)');
%                 ylabel('$/kWh');
%                 %axis([500,1500,0,0.02]);
%
%                 figure(2)
%                 plot(DH,TDamkW);
%                 title('Dam cost');
%                 xlabel('Dam height (m)');
%                 ylabel('$/kW');
%
%
%                 figure(3)
%                 plot(DH,DL);
%                 title('Dam height and dam width');
%                 xlabel('Dam height (m)');
%                 ylabel('Dam width (m)');
%                 %%
% h=figure(3);
% bar(1:20,[COETurR' COEDamR' COEOMR' COEPenR' COEElecR' COELake'],0.5,'stack');
% legend('Turbine', 'Dam', 'OM','Penstock','Electro','Lake')
% title('COE stacked');
% ylabel('$/kWh');
% xlabel('Dam width (*90m)');
%
% pathname = fileparts('Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\NAM\output\RD\');
% figfile = fullfile(pathname, sprintf('RDam_costperkWh_b%d_outlet%d.png',nbasin,outlet));
% saveas(h,figfile);
%
% %%
% hh=figure(4);
% bar(1:20,[AnTurR' AnDamR' AnOMR' AnPenR' AnElecR' AnLakeR'],0.5,'stack');
% legend('Turbine', 'Dam', 'OM','Penstock','Electro','Lake')
% title('Investments stacked');
% ylabel('Annualized investments');
% xlabel('Dam width *90m');
%
% pathname = fileparts('Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\NAM\output\RD\');
% figfile = fullfile(pathname, sprintf('RDam_Invest_b%d_outlet%d.png',nbasin,outlet));
% saveas(hh,figfile);

% end
%% Formula to re-calculate starting with SpecCapCost

% (((OptInv * AnnFac)/(P(b)*1e-3))  / (LF*8760*eta)) + TransCost
