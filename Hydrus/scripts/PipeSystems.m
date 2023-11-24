%% Pipe-systems procedure

%% Find pipe inlets
disp('Find Pinlets');

find_Pinlets

%% Determine lat lon coordinates of outlets and inlets
disp('Calculate pipelenght and lat lon');

Ppipelength

%% Running cost model
disp('Running pipe cost model');

Pcostmodel

%% Q-decreaser
disp('Q decreaser for pipe systems');

QPdecreaser

%% Store info per Q decreaser loop
disp('Store info per Q decreaser loop')

PCOEend{ndl}=aCOEPmin; % $/kWh Dam-Pipe COE
PPnetend{ndl}=aPPnetmin; % GWh Dam-Pipe energy potential
Platminend{ndl} = lat_Pin_min;
Plonminend{ndl}= lon_Pin_min;
rPinMinend{ndl} = rPinMin;
cPinMinend{ndl}= cPinMin;
aPInletminEnd{ndl}=aPInletMin;
aPinlet_windowMinEnd{ndl}=aPinlet_windowMin;
dfQPminEnd{ndl}=dfQPmin;
dfPLminEnd{ndl}=dfPLmin;
CostElementsPMinEnd{ndl}=CostElementsPMin;
Pinlet_winend{ndl}=Pinlet_win;
HeadPminend{ndl}=HeadPmin;
HeadraceLPminend{ndl}=HeadraceLPmin;
QDesignPinletMinend{ndl}=QDesignPinletMin;
QDesignLFPinletMinend{ndl}=QDesignLFPinletMin;
QDesignMeanPinletMinend{ndl}=QDesignMeanPinletMin;
dfPLminend{ndl}=dfPLmin;
ZPinletMinend{ndl} = ZPinletMin;
nPipePMinend{ndl} = nPipePMin;
OptInvPMinend{ndl} = OptInvPMin;
PPMinend{ndl} = PPMin;
OptSpecCapPMinend{ndl} = OptSpecCapPMin;
PSurfaceLake15sMinend{ndl} = PSurfaceLake15s;
accPminend{ndl} = accPmin;
OpthfPminend{ndl} = OpthfPMin;
OptDminend{ndl} = OptDPMin; %Optimal pipe diameter

%%
aDPPnetEnd4  = horzcat(PPnetend{:})'; % GWh  Dam-Pipe systems Pnet
fprintf('Total P potential: %0.0f GWh\n',sum(aDPPnetEnd4(~isnan(aDPPnetEnd4))));
