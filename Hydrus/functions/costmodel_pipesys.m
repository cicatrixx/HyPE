function [COEminRoR, TheoryPnet_GWh_min, Actual_PnetRoR_kWh_min, CostElements, npipemin, OptInv, OptP_W, OptSpecCap, Opthf, OptDP] = costmodel_pipesys(Zout,Zin,L,Qdesign,Qdesign_mean,LF,Dis2Trans,cost_constr,cost_lim,Hazardrate_RD_prct,Dis2Road,enablesmallDPcost,Dis2Settlement,LandValCost, ResettleCost, makechanges2cost)
% For each inlet-outlet with provided penstock length, compare the unit production
% cost for combination of 1-4 pipe of equal diameters rangeing from 0.1-4m 
% (10 possibilities tried out). Design flow is split across the pipes. 
% Use 10 iterations to calculate friction loss (hf) for all diameter and number of pipe combinations.
% Deselect diameters that have hf > available gross head. For valid diameters, 
% evaluate costs. Then pick cheapest combination of number of pipes and diameter. 

% ISSUE TO FIX: OptDP current returns idx for diameter rather than diameter
% value itself, need to fix that
run('config_costmodels.m')

% Make changes to config vals if any
if exist('makechanges2cost','var')
    %disp("Changes made to costconfig")
    eval(makechanges2cost);
end
%% Loop through npipe and pipe diameters
Hgross      = Zin-Zout;        %Gross head (m)
npipe_DP       = linspace(1,4,4);  % Number of pipes
D_DP           = linspace(0.1,4,10);  %Pipe diameter (meters)

for j=1:numel(npipe_DP)
    %% Physical calculations
    Atunnel = pi*(D_DP/2).^2;   %Cross sectional Tunnel(m2)
    f = 0.01;                %Initial friction number
    %Q = Qdesign;
    u = (Qdesign/npipe_DP(j))./Atunnel; % Water velocity (m/s)
    
    % Calculate friction factor in iterations
    for c = 1:10
        hf = f*L.*u.^2 ./ (2*g*D_DP);                            % Head loss
        Re = D_DP.*rho.*u./miu;                                   % Reynold number
        RHS = -2.*log10(e./(3.71.*D_DP) + 2.51./(Re.*sqrt(f)));  % RHS of Colebrook-White eq
        f = (1./RHS).^2;                                      % New friction factor
    end
    for i=1:numel(D_DP)
        if hf(i)>Hgross; hf(i)=NaN; end
    end
    
    % evaluate theory potential
    PTgross_MW = rho*g*Qdesign_mean*Hgross*1e-6;    % Theoretical gross potential in MW
    PTgross_GWh = rho*g*Qdesign_mean*Hgross*8760*1e-9;    % Theoretical gross potential
    
    %% Select tech/econ params for small or large project
    if enablesmallDPcost && PTgross_MW <50
        eta_generation_DP =eta_generation_DP_small;
        eta_transmission_DP=eta_transmission_DP_small;
        lifetime_DP=lifetime_DP_small;
        Ownersrate_DP  = Ownersrate_DP_small;
        OMshare_DP =OMshare_DP_small;
        DisPowerline_DP=Dis2Settlement;
        Hazardrate_DP=Hazard_DP_small;
        
        LandValCost=0; ResettleCost =0; % set to 0 for small DPs
    else
        eta_generation_DP =eta_generation_DP_large;
        eta_transmission_DP=eta_transmission_DP_large;
        lifetime_DP=lifetime_DP_large;
        Ownersrate_DP  = Ownersrate_DP_large;
        OMshare_DP =OMshare_DP_large;
        DisPowerline_DP=Dis2Trans;
        Hazardrate_DP=Hazard_DP_large_add+Hazardrate_RD_prct/100;
    end
    
    %% evaluate actual potential for multiple diameter options
    TheoryPnet_GWh = rho*g*Qdesign_mean*(Hgross-hf)*8760*1e-9; % Theoretical net potential
    Pror_W = eta_generation_DP*eta_transmission_DP*rho*g*(Hgross-hf).*Qdesign;            % Capacity based on design Q (W)
    
    Pror_W_end{j} = Pror_W;
    hf_end{j} = hf/Hgross;                            % Percentage loss due to head loss
    PnetRoR_kWh = (Pror_W*1e-3*8760).*LF;                   % Yearly production with LF (kWh)
    
    %% Cost calculations - components that are not optional
    %Different cost metrics based on capacity, Q, head and dam height
    %CostGunduz = -0.0254 + 1.81*QQ + 0.0848*(Hgross-Ahf) + 0.0553*L;    %Cost based on Gunduz & Sahin No transmission and 100y flood Q
    %CostORNL = 1.2 * 110168 * DH^-0.35*(P*1e-3).^-0.3;   %Initial capital cost ORNL formula ($/kW)
    TurIrena{j} = 1e6 *(1.1943* (Pror_W*1e-6).^0.7634) * IR2005;  %Turbine investment costs ($) IRENA hydropower, empirical R2=0.94
    %     Pelton2jV = (1105.2*QQ.^-0.5106).*(P*1e-3)*ER;        %Pelton turbine Veileder p156 2-jetter high head Q<10 h=1000
    %     Pelton6jV = (1559.6*QQ.^-0.5179).*(P*1e-3)*ER;        %Pelton turbine Veileder p156 6-jetter high head 10>Q<35 h=1000 ($)
    %     FrancisV = (1439.1*QQ.^-0.3143).*(P*1e-3)*ER;         %Francis turbine Veileder p157 0>Q<160 h=300 ($)
    %     KaplanV = (11730.5*QQ.^-0.2953).*(P*1e-3)*ER;         %Kaplan turbine Veileder p158 0>Q<400 h=15 ($)
    %     SpecTurcost = TurIrena./(P*1e-3);                     %Specific investment costs ($)
    
    %Total electrical mechanical costs based on Veileder p145
    PEM{j} = (3.9142*(Pror_W*1e-6).^0.6622) * ER_noktoUSD2010 *1e6;       %Total electro-technical equipment $
    PEMkW = PEM{j}./(Pror_W*1e-3);                          %Total electro-technical equipment $/kW
    
    %Steel penstock in tunnel costs based on Veileder report p94
    PenstockCost = (6*D_DP + 9.4)*1e3*ER_noktoUSD2010*Hgross; % $ Surface Penstock ($)
    
    %Powerstation surface Veileder p106
    PSSnok = -0.0006*Qdesign.^2 + 0.67*Qdesign - 6.95;    % mil NOK/m (head 10-40m, 1 powerunit)
    PSSnok(PSSnok<20) = 20;                     % Not lower than 20mil NOK/m
    PSS(j) = PSSnok*ER_noktoUSD2010*1e6;                        % $
    
    %Distribution Cost
    DisCostNOK = Powerline_allocator(DisPowerline_DP,Pror_W*1e-6); % NOK
    DisCost{j}=DisCostNOK*ER_noktoUSD2010; % $2010
    
    %% Components only for small or large DP
    if enablesmallDPcost && PTgross_MW <50
        % components only for small DP
        
        % Overground penstock
        TPCostL{j} = ( PenstockCost) * npipe_DP(j);
        
        % Small Desanding basin Intake Veileder 3(pg17)
        IntakeStructureNOK= -0.0033*Qdesign^2 + 0.1185*Qdesign + 0.3438;  % mil NOK/m
        Intake(j)= IntakeStructureNOK*ER_noktoUSD2010; % $2010
        
        RoadCost{j} = 0*TPCostL{j};
        Misc{j} =  0*TPCostL{j};
    else
        % components only for large DP
        %Tunnel blasted costs based on Veileder report p67
        CorrMP = 0.0054*(L*1e-3)^2 - 0.0039*(L*1e-3) + 0.9671; %lenght multiplier correlation
        PCost = (219.99*Atunnel + 13658) * ER_noktoUSD2010;   %Total price in ($/m2)
        PCostL = PCost * CorrMP;           %Cost times lenght multiplier ($/m)
        TPCostTunnel = PCostL*(L-Hgross);          %Total costs ($)
        TPCostL{j} = (TPCostTunnel + PenstockCost) * npipe_DP(j);
        
        %Underground powerstation added cost of blasted volume Veileder2 (p101) or Veileder3 (p18)
        BlastedV_m2= 78*Hgross^0.5*Qdesign^0.7*npowerunits^0.1;
        BlastCost=BlastedV_m2*2250*ER_noktoUSD2010;  % $
        PSS(j) = PSS(j)+BlastCost;                        % $
        PSSkW = PSS(j)./(Pror_W*1e-3);                   % $/kW
        
        % large Desanding basin Intake Veileder 3(pg17)
        IntakeStructureNOK = -0.0002*Qdesign^2 + 0.392*Qdesign + 1.1666;  % mil NOK/m
        Intake(j)= IntakeStructureNOK*ER_noktoUSD2010; % $2010
        
        %TransmissionCost
        TransCost=COETransCost*PnetRoR_kWh;
        DisCost{j}=DisCost{j}+TransCost;
        
        %RoadCost
        RoadCostNOK = road_allocator(Dis2Road,Pror_W*1e-6); % NOK
        RoadCost{j} =RoadCostNOK*ER_noktoUSD2010; % $2010
        
        % Fish pasage mitigation cost DOEwater p32: title: Estimation of Economic Parameters of U.S. Hydropower Resources
        FishCost = 1.3e6 *(Pror_W*1e-6).^0.56 * IR2002;
        
        %Miscellaneous
        Misc{j} = (-38.795.*log(Qdesign) + 309.89).*(Pror_W*1e-3)*ER_noktoUSD2010; %$
        
        % Add fish cost and provincial water use to misc costs
        Misc{j} = Misc{j} + FishCost + PWC_pak * ER_paktoUSD * PnetRoR_kWh;         
    end
    
    %% Calculate annulized costs
    AnnFac = interest/(1-((1+interest)^(-lifetime_DP)))  ; % Annuity factor using interest (10%) and economic lifetime (40 years) = 10,23%
    AnnualCost = AnnFac * (TurIrena{j}+ TPCostL{j} + PEM{j} + PSS(j)  + Intake(j)+ DisCost{j} + RoadCost{j} +  Misc{j} + LandValCost  + ResettleCost);
    COE = AnnualCost./PnetRoR_kWh;                           % Cost of electricity ($/kWh)
    
    % Collecting cost information
    AnTur = AnnFac * TurIrena{j};    %Turbine (annualized $)
    AnPen = AnnFac * TPCostL{j};     %Pipe (annualized $)
    AnElec = AnnFac * PEM{j};         %Electro-technical equipment ($)
    AnPS = AnnFac * PSS(j);           %Powerstation($)
    AnIn= AnnFac * Intake(j);           %Intake($)
    AnDis = AnnFac * DisCost{j};       %DisCost ($)
    AnRoad = AnnFac * RoadCost{j};       %RoadCost ($)
    AnMisc = AnnFac * Misc{j};        %Misc ($)
    AnLV = AnnFac * LandValCost;        %Misc ($)
    AnRC = AnnFac * ResettleCost;        %Misc ($)
    
    %
    AnOM = AnnualCost * OMshare_DP;         %OM cost as fraction of total investment costs
    AnCtot = AnTur + AnPen + AnElec + AnPS+ AnIn+ AnDis  + AnRoad +  AnMisc + AnOM + AnLV + AnRC;
    AnOwner = AnCtot * Ownersrate_DP;
    AnHazard = AnCtot * Hazardrate_DP;
    AnCtot2 = AnCtot + AnOwner + AnHazard;
    AnCtot_end{j} = AnCtot2;
    
    COETur = AnTur./PnetRoR_kWh;
    COEPen = AnPen./PnetRoR_kWh;
    COEElec = AnElec./PnetRoR_kWh;
    COEPS = AnPS./PnetRoR_kWh;
    COEIn = AnIn./PnetRoR_kWh;
    COEDis = AnDis./PnetRoR_kWh;
    COERoad = AnRoad./PnetRoR_kWh;
    COEMisc = AnMisc./PnetRoR_kWh;
    COELV = AnLV./PnetRoR_kWh;
    COERC = AnRC./PnetRoR_kWh;
    COEOM = AnOM./PnetRoR_kWh;
    COEOwner = AnOwner./PnetRoR_kWh;
    COEHazard = AnHazard./PnetRoR_kWh;
    
    COETot_np{j} = COETur + COEPen + COEElec + COEPS + COEIn + COEDis + COERoad + COEMisc + COELV + COERC + COEOM + COEOwner + COEHazard;
    
    %% Nan 0s for min search
    for jj=1:numel(COETot_np{j})
        if COETot_np{j}(jj)==0; COETot_np{j}(jj)=NaN;end
    end
    
    %% Find optimal pipe diameter for each npipe
    [COEminRoR_np(j) idx(j)]=min(COETot_np{j});
    
    Pror_W_min(j) = Pror_W_end{j}(idx(j));
    hfmin(j) = hf_end{j}(idx(j));
    AnCtotmin(j) = AnCtot_end{j}(idx(j));
    
    TheoryPnet_GWh_min_np(j) = TheoryPnet_GWh(idx(j));
    PnetRoR_kWh_min_np(j) = PnetRoR_kWh(idx(j));
    %
    CostElements_np{j}(1)=COETur(idx(j));
    CostElements_np{j}(2)=COEOM(idx(j));
    CostElements_np{j}(3)=COEPen(idx(j));
    CostElements_np{j}(4)=COEElec(idx(j));
    CostElements_np{j}(5)=COEPS(idx(j));
    CostElements_np{j}(6)=COEIn(idx(j));
    CostElements_np{j}(7)=COEMisc(idx(j));
    CostElements_np{j}(8)=COEDis(idx(j));
    CostElements_np{j}(9)=COEOwner(idx(j));
    CostElements_np{j}(10)=COEHazard(idx(j));
    CostElements_np{j}(11)=COELV(idx(j));
    CostElements_np{j}(12)=COERC(idx(j));
    
    TurIrenaMin(j) = TurIrena{j}(idx(j));
    TPCostLMin(j) = TPCostL{j}(idx(j));
    PEMMin(j) = PEM{j}(idx(j));
    PSSMin(j) = PSS(j);
    InMin(j) = Intake(j);
    DisCostMin(j) = DisCost{j}(idx(j));
    RoadCostMin(j) = RoadCost{j}(idx(j));
    MiscMin(j) = Misc{j}(idx(j));
    LVMin(j) = LandValCost;
    RCMin(j) = ResettleCost;
    
    if cost_constr==1
        if COEminRoR_np(j)>cost_lim; %Remove damsystem higher than x$/kWh
            COEminRoR_np(j)=NaN;
            PnetRoR_kWh_min_np(j)=NaN;
            TheoryPnet_GWh_min_np(j)=NaN;
            CostElements_np{j}=NaN;
            
        end
        %fprintf('no #%d, h #%d, in #%d: COE = %.2f $/kWh and Capacity = %.2f MW\n',k, heights, inletsp, COETotMx2, PMx2*1e-6)
        
    else
        %fprintf('no #%d, h #%d, in #%d: COE = %.2f $/kWh and Capacity = %.2f MW\n',k, heights, inletsp, COETotMx2, PMx2*1e-6)
        
    end
    %     fprintf('Power is %.2f MWh.\n',PnetRoRmin_np(j)*1e-6);
    %     fprintf('Lowest price is %.2f $/kWh.\n',COEminRoR_np(j));
    %     fprintf('Optimal pipe diameter %.2f m.\n\n',D(idx(j)));
    
    %% figure
    %     figure(10);
    %     plot(D,COETot_np{j},'color',cc(j,:), 'DisplayName',['Npipe = ',num2str(npipe(j))],'LineWidth',1.2);
    %     legend('-DynamicLegend')
    %     title('Production price');
    %     ylabel('$/kWh');
    %     xlabel('Pipe diameter');
    %     hold on
    %
    %     figure(2);
    %     bar(D,[COETur',COEOM', COEPipe',COEElec',COEPS', COEMisc', COEDis', COEOwner'],0.5,'stack')
    %     legend('Turbine','OM','Pipe','Electro','PS','Misc','Dis','Owner')
    %     legend('-DynamicLegend','location','Best')
    %     title('COE stacked','FontWeight','bold');
    %     ylabel('$/kWh');
    %     xlabel('Pipe diameter (m)');
    %
    %     %%
    %     figure(3);
    %     plot(D,PnetRoR*1e-6,'color',cc(j,:), 'DisplayName',['Npipe = ',num2str(npipe(j))],'LineWidth',1.2);
    %     legend('-DynamicLegend')
    %     title('Energy production');
    %     ylabel('MWh');
    %     xlabel('Pipe diameter');
    %     hold on
    
end

%% Find optimal npipes
[COEminRoR ind] = min(COEminRoR_np); % $/kWh
CostElements = CostElements_np{ind};
TheoryPnet_GWh_min = TheoryPnet_GWh_min_np(ind); %theory GWh
Actual_PnetRoR_kWh_min = PnetRoR_kWh_min_np(ind); %actual kWh
npipemin = npipe_DP(ind);
OptP_W = Pror_W_min(ind); %W
OptInv = AnCtotmin(ind); % $ annualized
Opthf = hfmin(ind); % percentage of gross head lost due to friction loss (head loss)
%OptSpecCap = OptInv/(OptP*1e-3); % $/kW Annualized speccapcost
OptDP = idx(ind);

OptInv1 = TurIrenaMin(ind) + TPCostLMin(ind) + PEMMin(ind) + PSSMin(ind) ...
    +InMin(j)+ DisCostMin(ind) +RoadCostMin (ind)+ MiscMin(ind) + LVMin(ind) + RCMin(ind); % $ absolute investments
OptInvOM = OMshare_DP * OptInv1;
OptInvsOwners = (OptInv1 + OptInvOM) * Ownersrate_DP;
OptInvsQuake = (OptInv1 + OptInvOM) * Hazardrate_DP;
OptInv2 = OptInv1 + OptInvOM + OptInvsOwners + OptInvsQuake;

OptSpecCap = OptInv2/(OptP_W*1e-3); % $/kW not annualized

