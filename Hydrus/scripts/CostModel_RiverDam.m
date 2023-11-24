%% Costmodel_RiverDam
% Calculate cost and potential of river dam

COETotRD = 0;
RDPnet=0;
clear DH
for  k = 1:no
    
    if isnan(ro(k))==1;
        COETotRD(k)=0;
        RDPnet(k)=0;
        OptDH(k)=0;
        OptDL(k)=0;
        RDCostElements{k}=0;
        continue; end; % Skip empty arrays
    
    if numel(dfldamRD{k})==0;
        COETotRD(k)=0;
        RDPnet(k)=0;
        OptDH(k)=0;
        OptDL(k)=0;
        RDCostElements{k}=0;
        continue; end; % Skip empty arrays
    
    if numel(dfhdamRD{k})==0;
        COETotRD(k)=0;
        RDPnet(k)=0;
        OptDH(k)=0;
        OptDL(k)=0;
        RDCostElements{k}=0;
        continue; end; % Skip empty arrays
    
    DH= single(dfhdamRD{k});
    
    if MiniHydro_special==1 & Qoutlets_design(k) < bighydro_cutoff
        [COETotRD(k) RDPnet(k) OptP(k) OptDH(k) OptDL(k) RDCostElements{k}]=MiniCostmodel(dfldamRD{k},DH,PopCost{k},Qoutlets_design(k),Qoutlets_design_LF(k),DisOutlet(k),LandValueRDlake{k},RDDepth(k),nbasin,k,cost_constr,cost_lim);
    else
        [COETotRD(k) RDPnet(k) OptP(k) OptDH(k) OptDL(k) RDCostElements{k} OptInv(k) RDP(k) OptPop(k) OptLV(k) OptSpecCap(k)]=RDcostmodel(dfldamRD{k},DH,PopCost{k},Qoutlets_design(k),Qoutlets_design_LF(k),DisOutlet(k),LandValueRDlake{k},RDDepth(k),nbasin,k,cost_constr,cost_lim,Quakerate(k));
    end
    
    PopDisplacedOpt(k) = PopDisplaced{k}(OptDH(k));
end


%% To remove zeros
for k =1:no
    if COETotRD(k)==0; COETotRD(k)=NaN; end;
    if RDPnet(k)==0; RDPnet(k)=NaN;end;
end

COETotRD = COETotRD';
RDPnet = RDPnet'; % kWh
RDPnet = RDPnet*1e-6; % GWh


