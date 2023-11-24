%% Find Pinlet locations

for k = 1:numel(outlets);
    
    if isnan(ro(k))==1; Pinlet_win{k}=[]; continue; end;
    
    %Data windows are created in find_outlet_wins
    
    [nr_ZiW,nc_ZiW] = size(Z_inlet_win{k});
    
    Z_inlet_win_temp=0;
    rr1=inlet_win+1-sradius;
    rr2=inlet_win+1+sradius;
    cc1=inlet_win+1-sradius;
    cc2=inlet_win+1+sradius;
    Z_inlet_win_temp = Z_inlet_win{k}(rr1:rr2,cc1:cc2);
    
    Zmin = Z_inlet_win{k}(inlet_win+1,inlet_win+1);
    Zmax = max(Z_inlet_win_temp(:));
    Zupdown = Zmax-Zmin;
    ZPin = linspace(Zmin+(0.1*Zupdown),Zmax-(0.1*Zupdown),ni);
    
    Zout=  sub2ind(size(Z_inlet_win{k}),inlet_win+1,inlet_win+1); %index of center dam outlet location
    
    %figure(2);clf;imagesc(Z_inlet_win{k});axis image; colormap(flipud(gray));hold on
    %plot(inlet_win+1,inlet_win+1,'.r','markersize',20); hold off
    
    for l = 1:numel(ZPin)
        
        %Pinletsw = inlet_finder(Z_inlet_win{k},Q_inlet_win{k},acc_inlet_win{k},fdir_inlet_win{k},WDPA_PL10_inlet_win{k},WDPA_PL20_inlet_win{k},protarea_constr,nr_ZiW,nc_ZiW,Zout,ZPin(l),deselect_mainstream,index_inlet_win,adir_inlet_win{k});
        Pinletsw = inlet_finder(Z_inlet_win{k},Q_inlet_win{k},acc_inlet_win{k},fdir_inlet_win{k},WDPA_PL10_inlet_win{k},WDPA_PL20_inlet_win{k},protarea_constr,nr_ZiW,nc_ZiW,Zout,ZPin(l),deselect_mainstream,index_inlet_win,adir_inlet_win{k},flowdist_inlet_win{k},sradius);
        
        Pinlet_win{k}{l} = find(Pinletsw);
        %convert inlet maps to arrays of rows/columns
        rPin_w{k}{l} = rem(Pinlet_win{k}{l}-1, nr_ZiW)+1;
        cPin_w{k}{l} = fix((Pinlet_win{k}{l}-1)/nr_ZiW)+1;
        % Converting r,c to Columbia catchment r,c
        rPin{k}{l} = ro(k)+(rPin_w{k}{l}-inlet_win-1);
        cPin{k}{l} = co(k)+(cPin_w{k}{l}-inlet_win-1);
        Pinlet{k}{l} = sub2ind(size(Z),rPin{k}{l},cPin{k}{l});
        ZPinlet{k}{l} = Z(Pinlet{k}{l});
        QPinlet{k}{l} = Q(Pinlet{k}{l});
        accPinlet{k}{l} = acc(Pinlet{k}{l});
        QDesignPinlet{k}{l} = Qdesign(Pinlet{k}{l});
        QDesignLFPinlet{k}{l} = Qdesign_LF(Pinlet{k}{l});
        QDesignMeanPinlet{k}{l} = Qdesign_mean(Pinlet{k}{l});
        PRegion_id{k}{l} = Regions(rPin{k}{l},cPin{k}{l});
        PCountry_id{k}{l} = Countries(rPin{k}{l},cPin{k}{l});

        
        if slackflow_constraint==1
            QDesignPinlet{k}{l} = QDesignPinlet{k}{l} *( 1-slackflow);
        end
        
        %Check
        %figure(3);clf;imagesc(log(Q_inlet_win{3}));axis image; colormap(flipud(gray));hold on
        %plot(inlet_win+1, inlet_win+1,'b.','markersize',24);
        %plot(cPin_w{3}{1}, rPin_w{3}{1},'r.','markersize',24);
        %hold off
        
        %figure(3);clf;imagesc(log(Q));axis image; colormap(flipud(gray));hold on
        %plot(cPin{3}{1}, rPin{3}{1},'r.','markersize',24);
        %hold off
        
    end
end

%% Variable creator


for k = 1:numel(outlets)
    
    COEP{k}=[];
    PL{k}=[];
    latP{k}=[];
    lonP{k}=[];
    PPnet{k}=[];
    PPipeDia{k}=[];
    OptInvP{k}=[];
    PP{k}=[];
    bMinElevP(k)=0;
    aCOEPmin(k)=0;
    aPPnetmin(k)=0;
    aCOEPmin(k)=NaN;
    aPPnetmin(k)=NaN;
    CostElementsPMin{k}=[];
    dfQPmin(k)=NaN;
    accPmin(k)=NaN;
    dfPLmin(k)=NaN;
    HeadPmin(k)=NaN;
    HeadraceLPmin(k)=NaN;
    aPInletMin(k)=NaN;
    aPinlet_windowMin(k)=NaN;
    lat_Pin_min(k)=NaN;
    lon_Pin_min(k)=NaN;
    rPinMin(k)=NaN;
    cPinMin(k)=NaN;
    QDesignPinletMin(k) =NaN;
    QDesignLFPinletMin(k) =NaN;
    QDesignMeanPinletMin(k) = NaN;
    ZPinletMin(k) =NaN;
    nPipePMin(k) =NaN;
    OptInvPMin(k) =NaN;
    PPMin(k) =NaN;
    DeselectedSites_unSortRD{k}=[];
    DeselectedGrandSites_unSortRD{k}=[];
    MiniFlag(k) = 0;
    OptInv(k) =NaN;
    OptSpecCap(k) =NaN;
    RDP(k) = NaN;
    RDVolumeLake15s(k) = NaN;
    RDSurfaceLake15s(k) = NaN;
    OptPop(k) = NaN;
    OptLV(k) = NaN;
    OptSpecCapP{k} = [];
    OpthfP{k} = [];
    OptDP{k} = [];
    OptSpecCapPMin(k) = NaN;
    OpthfPMin(k) = NaN;
    PSurfaceLake15s(k) = 1e-3;
    OptDPMin(k) = NaN;
    
    if numel(Pinlet_win{k})==0; continue; end;
    
    for l = 1:ni
        COEP{k}{l}=[];
        PL{k}{l}=[];
        latP{k}{l}=[];
        lonP{k}{l}=[];
        PPnet{k}{l}=[];
        PPipeDia{k}{l}=[];
        OptInvP{k}{l}=[];
        PP{k}{l}=[];
        COEPminElev{k}(l)=0;
        bMinP{k}(l)=0;
        OptSpecCapP{k}{l}=[];
        OpthfP{k}{l}=[];
        OptDP{k}{l}=[];
        
        for m = 1:numel(Pinlet_win{k}{l})
            
            COEP{k}{l}(m)=NaN;
            PL{k}{l}(m)=NaN;
            latP{k}{l}(m)=NaN;
            lonP{k}{l}(m)=NaN;
            PPnet{k}{l}(m)=NaN;
            PPipeDia{k}{l}(m)=NaN;
            OptInvP{k}{l}(m)=NaN;
            PP{k}{l}(m)=NaN;
            OptSpecCapP{k}{l}(m)=NaN;
            OpthfP{k}{l}(m)=NaN;
            OptDP{k}{l}(m)=NaN;
        end
    end
end

for k=1:nd
    for l = 1:numel(outlets)
        for m = 1:ni
            DeselectedSites_unSortDP{k}{l}{m}=[];
        end
    end
end

DeselectedSites = 0;
DeselectedSites_final = {};
