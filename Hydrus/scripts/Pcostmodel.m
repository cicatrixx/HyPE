%% Pipe cost model

for  k = 1:no
    if isnan(ro(k))==1; continue; end;
    for l = 1:ni
        for m = 1:numel(Pinlet_win{k}{l})
            
            if isnan(lonP{k}{l}(m))==1;
                COEP{k}{l}(m)= NaN;
                PPnet{k}{l}(m) = NaN;
                CostElementsP{k}{l}{m} = NaN;
                nPipeP{k}{l}(m)=NaN;
                continue; end; 
            
            [COEP{k}{l}(m) PPnet{k}{l}(m) CostElementsP{k}{l}{m} nPipeP{k}{l}(m) OptInvP{k}{l}(m) PP{k}{l}(m) OptSpecCapP{k}{l}(m) OpthfP{k}{l}(m) OptDP{k}{l}(m)] = costmodel_pipesys(Zoutlets(k),ZPinlet{k}{l}(m),PL{k}{l}(m),QDesignPinlet{k}{l}(m),QDesignMeanPinlet{k}{l}(m),QDesignLFPinlet{k}{l}(m),DisOutlet(k),cost_constr,cost_lim,Quakerate(k));
            
            if COEP{k}{l}(m)==0; COEP{k}{l}(m)=NaN;end;
            if PPnet{k}{l}(m)==0; PPnet{k}{l}(m)=NaN;end;
            PPnet{k}{l}(m) = PPnet{k}{l}(m) *1e-6; % kWh to GWh
            
        end 
    end
end

%% Select Cheapest
% First find lowest COE per elevation level
for k=1:no
    
    if isnan(ro(k))==1; continue; end;
    
    for l=1:ni
        if isempty(COEP{k}{l})==1; continue; end;
        [COEPminElev{k}(l),bMinP{k}(l)] = min(COEP{k}{l}); %Lowest index per elevation level
    end
end

for k=1:no
     
    if isnan(ro(k))==1; continue; end;
     
    bMinP{k}(find(bMinP{k}==0))=NaN;
    COEPminElev{k}(find(COEPminElev{k}==0))=NaN;
end

% Secondly find lowest COE among elevation levels
for k=1:no
    
    if isnan(ro(k))==1; continue; end;
    if isempty(COEPminElev{k})==1; continue; end;
    
    
    [~,bMinElevP(k)] = min(COEPminElev{k}); %Lowest index of all elevation levels
end

%% Transfer min info into actual data
clear a1 a2
for k= 1:no
    
    if isnan(ro(k))==1; continue; end;
    if bMinElevP(k)==0; continue; end;
    if isempty(bMinP{k})==1; continue; end;
    
    a1(k)=bMinElevP(k); % Elevation level
    a2(k)=bMinP{k}(bMinElevP(k)); %inlet min
    
    if a2(k)==0; continue; end;
    if isnan(a2(k))==1; continue; end;
    
    aCOEPmin(k) = COEP{k}{a1(k)}(a2(k));
    aPPnetmin(k) = PPnet{k}{a1(k)}(a2(k));
    nPipeminP(k) = nPipeP{k}{a1(k)}(a2(k));
    
    
    if isempty(CostElementsP{k}{a1(k)})==1;
        CostElementsPMin{k} = 0;
        dfQPmin(k)= 0;
        dfPLmin(k)= 0;
        HeadPmin(k) = 0;
        aPInletMin(k) = 0;
        aPinlet_windowMin(k) = 0;
        lat_Pin_min(k) = 0;
        lon_Pin_min(k) = 0;
        rPinMin(k) = 0;    
        cPinMin(k) = 0;
        continue; end

    CostElementsPMin{k} = CostElementsP{k}{a1(k)}{a2(k)};
    dfQPmin(k)= QPinlet{k}{a1(k)}(a2(k));
    dfPLmin(k)= PL{k}{a1(k)}(a2(k));
    HeadPmin(k) = ZPinlet{k}{a1(k)}(a2(k)) - Zoutlets(k);
    HeadraceLPmin(k)= PL{k}{a1(k)}(a2(k)) - (ZPinlet{k}{a1(k)}(a2(k)) - Zoutlets(k));
    aPInletMin(k) = Pinlet{k}{a1(k)}(a2(k));
    aPinlet_windowMin(k) = Pinlet_win{k}{a1(k)}(a2(k));
    lat_Pin_min(k) = latP{k}{a1(k)}(a2(k));
    lon_Pin_min(k) = lonP{k}{a1(k)}(a2(k));
    rPinMin(k) = rem(aPInletMin(k)-1, nrw)+1;  %best inlet row
    cPinMin(k) = fix((aPInletMin(k)-1)/nrw)+1;  %best inlet col
    QDesignPinletMin(k) = QDesignPinlet{k}{a1(k)}(a2(k));
    QDesignLFPinletMin(k) = QDesignLFPinlet{k}{a1(k)}(a2(k));
    QDesignMeanPinletMin(k) = QDesignMeanPinlet{k}{a1(k)}(a2(k));
    ZPinletMin(k) = ZPinlet{k}{a1(k)}(a2(k));
    nPipePMin(k) = nPipeP{k}{a1(k)}(a2(k));
    OptInvPMin(k) = OptInvP{k}{a1(k)}(a2(k));
    PPMin(k) = PP{k}{a1(k)}(a2(k));
    OptSpecCapPMin(k) = OptSpecCapP{k}{a1(k)}(a2(k));
    accPmin(k) = accPinlet{k}{a1(k)}(a2(k));
    OpthfPMin(k) = OpthfP{k}{a1(k)}(a2(k));
    OptDPMin(k) = OptDP{k}{a1(k)}(a2(k));
   
end
