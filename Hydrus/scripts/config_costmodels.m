%% Generic
g           = 9.8;             %Gravitational acceleration (m/s2)
rho         = 1000;            %Density of water (kg/m3)
npowerunits=2;                  % number of power units in power station for large DP and R projects for calculating underground power station costs

%Generic cost factors
PWC_pak = 0.425;                %Provincial water use charge as per GoP 1984 (Pak/kWh)
interest    = 10/100;             %Interest on capital - discount rate
ER_paktoUSD = 0.0065 ;          %Exchange rate PAK to US$ 21/04/2021
ER_noktoUSD2010 = 0.1725;          %Exchange rate NOK to US$ 01/01/2010
IR2002      = 1.24;            %Inflation rate conversion t0 2010 USD: World Bank Real effective exchange rate index (2010 = 100)
IR2005      = 1.09;            %Inflation rate conversion t0 2010 USD: http://data.worldbank.org/indicator/PX.REX.REER?page=2
IR2013      = 0.99;            %Inflation rate conversion to 2010 USD: http://data.worldbank.org/indicator/PX.REX.REER
COETransCost   = 0.0034 * IR2013; %Fixed transmission cost per delivered kWh: http://www.eia.gov/forecasts/aeo/tables_ref.cfm Table 8 EIA ANNUAL ENERGY OUTLOOK 2015

%% DP costmodel constants
miu          = 0.001;           %Fluid viscosity of water (N-s/m2))
e           = 0.2;             %Roughness constant rough concrete (m)

%DP-large-specific
eta_generation_DP_large = 85/100;          %Efficiency turbine
eta_transmission_DP_large = 85/100;   %Transmission and distribution efficiency reported for NEPRA 2019
lifetime_DP_large    = 30;              %Years
Ownersrate_DP_large  = 10/100;            %Owners cost due to lead times IRENA 2012 Hydropower
OMshare_DP_large     = 3/100;          %As share of total investments
Hazard_DP_large_add     = 1/100;          %Added cost for large DP added on top of the RP rate as share of total investments

%DP-small-specific
eta_generation_DP_small = 80/100;          %Efficiency turbine for Madakhalast SHP
eta_transmission_DP_small  = 96/100;   %Transmission and distribution efficiency for Madakhalast SHP
lifetime_DP_small     = 15;              %Years
Ownersrate_DP_small   = 0;            %Owners cost due to lead times IRENA 2012 Hydropower
OMshare_DP_small      = 3.5/100;          %As share of total investments
Hazard_DP_small     = 0;          %As share of total investments

%% RDcostmodel constants
D_RD          = 3;              %Pipe diameter (meters)
eta_generation_RD = 90/100;     %Water to wire efficiency (turbine losses, pipe friction losses)
eta_transmission_RD = 85/100;   %Transmission and distribution efficiency reported for NEPRA 2019
lifetime_RD    = 40;            %Years
Ownersrate_RD  = 20/100;        %Owners cost due to lead times IRENA 2012 Hydropower
OMshare_RD     = 2.5/100;        %As share of total investments
%%disp("Cost parameters loaded")
