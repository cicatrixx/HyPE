%% Dam selector
% Routine to de-selects dams according to a certain priority

if sum(isnan(COETotRD))==numel(ro)
    RDlake_Opt=0;
    CheapestDam=0;
    DeselectedSites=0;
    DeselectedSites_unSort=0;
    COETotRDs=0;
    RDPnets=0;
end

%% Collecting all dam locations DP, P and RD

if runDP==1
    %Zeros are NaN
    aInletminEnd{1}(aInletminEnd{1}(:)==0)=NaN;
    aInletminEnd{2}(aInletminEnd{2}(:)==0)=NaN;
    outlets(outlets(:)==0)=NaN;
    
    DPlocs = horzcat(aInletminEnd{:})';
    Plocs  = horzcat(aPInletminEnd{:})';
    RDlocs = outlets';
    
    Damlocs = [DPlocs;Plocs;RDlocs]; %In terms of basin indices
    
    % Vector to keep track of type of systems
    SysID  = zeros((numel(RDlocs)+numel(DPlocs)+numel(Plocs)),1);
    SysID(1:numel(DPlocs))=1;                                      %DP systems
    SysID((numel(DPlocs)+1):(numel(DPlocs)+numel(Plocs)))=2;       %P systems
    SysID((numel(DPlocs)+numel(Plocs)+1):end)=3;                   %RD systems
    
    % Rows and columns
    ros = [horzcat(rinMinend{:})'; horzcat(rPinMinend{:})'; ro];
    cos = [horzcat(cinMinend{:})'; horzcat(cPinMinend{:})'; co];
    
    lats = [horzcat(latminend{:})'; horzcat(Platminend{:})'; latOut];
    lons = [horzcat(latminend{:})'; horzcat(Plonminend{:})'; lonOut];
    
    % %%
    % figure(1);clf;imagesc(log(Q));colormap(flipud(gray));axis image;
    % hold on
    % plot(co(find(~isnan(COETotRD))),ro(find(~isnan(COETotRD))),'.r','markersize',20)
    % plot(cinMinend{1}(find(~isnan(COEend{1}))),rinMinend{1}(find(~isnan(COEend{1}))),'.b','markersize',20)
    % plot(cPinMinend{1}(find(~isnan(PCOEend{1}))),rPinMinend{1}(find(~isnan(PCOEend{1}))),'.g','markersize',20)
    % hold off
    
else
    %Zeros are NaN
    outlets(outlets(:)==0)=NaN;
    
    Plocs  = horzcat(aPInletminEnd{:})';
    RDlocs = outlets';
    
    Damlocs = [Plocs;RDlocs]; %In terms of basin indices
    Grandlocs = GrandIdx; %Existing locations
    
    % Vector to keep track of type of systems
    SysID  = zeros((numel(RDlocs)+numel(Plocs)),1);
    SysID((1:numel(Plocs)))=1;       %P systems
    SysID((numel(Plocs)+1):end)=2;   %RD systems
    
    % Rows and columns
    ros = [horzcat(rPinMinend{:})'; ro];
    cos = [horzcat(cPinMinend{:})'; co];
    
    lats = [horzcat(Platminend{:})'; latOut'];
    lons = [horzcat(Plonminend{:})'; lonOut'];
    
end

%% Show
% [rRD,cRD]=ind2sub(size(acc),RDlocs);
% % [rDP,cDP]=ind2sub(size(acc),DPlocs);
% [rP,cP]=ind2sub(size(acc),Plocs);
% 
% figure(1);clf;imagesc(log(Q));axis image;colormap(flipud(gray));
% hold on
% % plot(cRD,rRD,'r.','markersize',10)
% % plot(cDP,rDP,'b.','markersize',20)
% % plot(cP,rP,'g.','markersize',10)
% plot(coss,ross,'g.','markersize',20)
% plot(c_damst,r_damst,'b.','markersize',20)
% hold off
%% Collecting COEs and Pnets

if runDP==1
    % COEs
    COEendDP = horzcat(COEend{:});
    COEendDP(COEendDP(:)==0)=NaN;
    
    COEendP = horzcat(PCOEend{:});
    COEendP(COEendP(:)==0)=NaN;
    
    COEAll = [COEendDP';COEendP';COETotRD];
    
    % Pnets
    PnetendDP = horzcat(DPPnetend{:});
    PnetendDP(PnetendDP(:)==0)=NaN;
    
    PnetendP = horzcat(PPnetend{:});
    PnetendP(PnetendP(:)==0)=NaN;
    
    PnetAll = [PnetendDP';PnetendP';RDPnet];
    
else
    % COEs
    COEendP = horzcat(PCOEend{:});
    COEendP(COEendP(:)==0)=NaN;
    
    COEAll = [COEendP';COETotRD];
    
    % Pnets
    PnetendP = horzcat(PPnetend{:});
    PnetendP(PnetendP(:)==0)=NaN;
    
    PnetAll = [PnetendP';RDPnet];
    
end

%% First, recalculate lake based on optimal dam height and check which lake floods which dams
%% DP-systems

if runDP==1
    if sum(isnan(COEendDP))~=numel(COEendDP)
        for ndd=1:2
            for k=1:no
                fprintf('Recalculating DP lakes outlet #%d of %d\n',k,no);
                
                clear DPupstream a OutletDeselect
                if isnan(ro(k))==1;
                    DeselectedSites_unSortDP{ndd}{k}=0;
                    continue; end;
                if isnan(aInletminEnd{ndd}(k))==1;
                    DeselectedSites_unSortDP{ndd}{k}=0;
                    continue; end;
                
                DPupstream = fastfindupstream_lim(acc,fdir,drow,dcol,aInletminEnd{ndd}(k));
                
                DPlake_Opt = DPupstream & Z < (Z(aInletminEnd{ndd}(k))-DPDepth{ndd}(k)+ahdammin{ndd}(k));
                
                LakeIdx = find(DPlake_Opt);
                a = ismember(Damlocs,LakeIdx);
                OutletDeselect = find(a);
                kc=k;
                if ndd==2; kc=k+no;end % Making sure it doesn't select itself in the second round
                DeselectedSites_unSortDP{ndd}{k} = OutletDeselect(OutletDeselect~=kc); % Which dams are flooded (unsorted)
                
                clear a
                a = ismember(Grandlocs,LakeIdx);
                GrandOutletDeselectDO = find(a);
                
            end
        end
    end
end

%% Show stuation before deselection

% for i=1:numel(DPlake_Opt{1})
%     a(i)=isempty(DPlake_Opt{1}{i});
% end
% b=find(~a);
%
% %
% figure(1);clf;
% img{1}= truecolorsc(DPlake_Opt{1}{b(1)},flipud(gray));
% %
% for i=1:numel(b)
%     if i==1; continue; end
%         img{i} = burnmask(img{i-1}, ~DPlake_Opt{1}{b(i)});
% end
% image(img{end}); axis image;
% %
% hold on
% plot(cinMinend{1}(:),rinMinend{1}(:),'.r','markersize',15)
% plot(cPinMinend{1}(:),rPinMinend{1}(:),'.g','markersize',15)
% plot(cinMinend{2}(:),rinMinend{2}(:),'.r','markersize',15)
% plot(cPinMinend{2}(:),rPinMinend{2}(:),'.g','markersize',15)
% plot(co,ro,'.b','markersize',15)
% hold off


%% RD-systems

if sum(isnan(COETotRD))~=numel(ro)
    for k=1:no
        fprintf('Recalculating lakes outlet #%d of %d\n',k,no);
        
        clear a OutletDeselect Z_upstream RDupstream dz_h RDlake_Opt
        if isnan(ro(k))==1;
            DeselectedSites_unSortRD{k}=0;
            DeselectedGrandSites_unSortRD{k}=0;
            continue; end;
        if isnan(COETotRD(k))==1;
            DeselectedSites_unSortRD{k}=0;
            DeselectedGrandSites_unSortRD{k}=0;
            continue; end;
        
        RDupstream = fastfindupstream_lim(acc,fdir,drow,dcol,outlets(k));
        
        RDlake_Opt = RDupstream & Z < (Zoutlets(k)-RDDepth(k)+OptDH(k));
        
        %For potential dams
        LakeIdx = find(RDlake_Opt);
        a = ismember(Damlocs,LakeIdx);
        OutletDeselect = find(a);
        if runDP==1
            DeselectedSites_unSortRD{k} = OutletDeselect(OutletDeselect~=(numel(DPlocs)+numel(Plocs)+k)); % Which dams are flooded (unsorted)
        else
            DeselectedSites_unSortRD{k} = OutletDeselect(OutletDeselect~=(numel(Plocs)+k)); % Which dams are flooded (unsorted)
        end
        
        %For existing dams
        clear a
        a = ismember(Grandlocs,LakeIdx);
        GrandOutletDeselectRD = find(a); %
        DeselectedGrandSites_unSortRD{k} = GrandOutletDeselectRD; % Which Grand dams are flooded (unsorted)
        
        %% Lake volume calculator m3
        Z_upstream = Z(RDupstream);
        dz_h = Zoutlets(k)-RDDepth(k)+OptDH(k)-Z_upstream;
        dz_h = max(0, dz_h);
        RDVolumeLake15s(k)  = sum(dz_h * 450 * (450*cosd(latOut(k)))); %Volume lake m3
        RDSurfaceLake15s(k) = max(1,sum(RDlake_Opt(:))) * 450 * (450*cosd(latOut(k))); %Surface Reservoir m2
    end
end

%% Collecting all deselected DamLocs
%clear DeselectedSites CheapestDam DeselectedSites_final

if runDP==1
    DeselectedSites_unSortDPP = horzcat(DeselectedSites_unSortDP{:});
    
    for k = 1:numel(DeselectedSites_unSortDPP); DeselectedSites_unSortP{k}=[]; end %Creating empty P array (floods nothing)
    
    DeselectedSites_unSort = [DeselectedSites_unSortDPP,DeselectedSites_unSortP,DeselectedSites_unSortRD];
    
else
    for k = 1:(nd*numel(DeselectedSites_unSortRD)); DeselectedSites_unSortP{k}=[]; end %Creating empty P array (floods nothing)
    
    DeselectedSites_unSort = [DeselectedSites_unSortP,DeselectedSites_unSortRD];
    DeselectedGrandSites_unSort = [DeselectedSites_unSortP,DeselectedGrandSites_unSortRD];
    
    % Lake surfaces
    LakeSurfacesAll = [horzcat(PSurfaceLake15sMinend{:})'; RDSurfaceLake15s'];
    
    % Flow accumulations
    accAll = [horzcat(accPminend{:})'; Accoutlets'];
    
    % Basin IDs
    BIDall = zeros(size(accAll));
    BIDall(:) = nbasin;
    
    % Continent IDs
    CIDall = zeros(size(accAll));
    CIDall(:) = CID;
end

%% Show stuation before deselection
% figure(1);clf;
%
% img{1}= truecolorsc(RDlake_Opt{1},flipud(gray));
%
% for i=1:numel(RDlake_Opt)
%     if i==1; continue; end
%
%     img{i} = burnmask(img{i-1}, ~RDlake_Opt{i});
% end
% image(img{end}); axis image;
%
% hold on
% plot(co(1:i),ro(1:i),'.r','markersize',15)
% hold off

%% Second, sorting through indexing

if sum(isnan(COETotRD))~=numel(ro)
    if exist('DeselectedSites_unSort')
        
        outletsRel=find(~isnan(COEAll)); %Find relevant outlet idx
        
        %Cost Priority
        minCOE = min(COEAll);
        COEIndexing = 1./(COEAll/minCOE); %Lowest COE gets 1
        if sum(isnan(COEIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in COE deselector');return;end
        if sum(isinf(COEIndexing(outletsRel))) > 0; returnID=1;disp('Infs in COE deselector');return;end
        
        %Power priority
        maxPnet = max(PnetAll);
        PnetIndexing = PnetAll/maxPnet; %Highest Pnet gets 1
        if sum(isnan(PnetIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in Pnet deselector');return;end
        if sum(isinf(PnetIndexing(outletsRel))) > 0; returnID=1;disp('Infs in Pnet deselector');return;end
        
        %Cost per power priority
        COEperPnet = COEAll./PnetAll;
        minCOEperPnet = min(COEperPnet);
        COEperPnetIndexing = 1./(COEperPnet/minCOEperPnet); %Lowest COEperPnet gets 1
        if sum(isnan(COEperPnetIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in COEPnet deselector');return;end
        if sum(isinf(COEperPnetIndexing(outletsRel))) > 0; returnID=1;disp('Infs in COEPnet deselector');return;end
        
        %Lake surface priority
        MaxSurfaceLake = max(LakeSurfacesAll);
        LakeIndexing = 1./(LakeSurfacesAll/MaxSurfaceLake); %Lowest surface gets 1
        if sum(isnan(LakeIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in Lake deselector');return;end
        if sum(isinf(LakeIndexing(outletsRel))) > 0; returnID=1;disp('Infs in Lake deselector');return;end
        
        %Flow accumulation priority
        MaxaccAll = max(accAll);
        AccIndexing = 1./(accAll/MaxaccAll); %Lowest flow acc gets 1
        if sum(isnan(AccIndexing(outletsRel))) > 0; returnID=1;disp('NaNs in acc deselector');return;end
        if sum(isinf(AccIndexing(outletsRel))) > 0; returnID=1;disp('Infs in acc deselector');return;end
        
        WeightedCOEPnet = (COEIndexing*betaC)+(PnetIndexing*betaP)+(COEperPnetIndexing*betaCP)+(LakeIndexing*betaL)+(AccIndexing*betaA);
        
        %Sorting
        [~, sortIdx] = sort(WeightedCOEPnet(outletsRel),'descend');
        
        CheapestDam = outletsRel(sortIdx);
        DeselectedSites = DeselectedSites_unSort(outletsRel(sortIdx)); %Sorted flooded dams
        
    else
        RDlake_Opt=0;
        CheapestDam=0;
        DeselectedSites=0;
        COEAll=0;
        PnetAll=0;
    end
end

%% Third, deselection starting with cheapest
if iscell(DeselectedSites)
    CheapestDam2=CheapestDam;
    for i=1:numel(CheapestDam2)
        if sum(ismember(CheapestDam2(i+1:end), DeselectedSites{i}))>=1;
            idxdeselect{i} = find(ismember(CheapestDam2(i+1:end), DeselectedSites{i})) + i;
            CheapestDam2(idxdeselect{i})=1;
            continue; end
    end
end

%% Fourth, deselection starting with most expensive
if iscell(DeselectedSites)
    CheapestDam3=CheapestDam2;
    for i=0:numel(CheapestDam3)-1
        if i==numel(CheapestDam3)-1; continue; end
        if sum(ismember(CheapestDam3(:), DeselectedSites{numel(CheapestDam3)-i})) >=1
            CheapestDam3(numel(CheapestDam3)-i)=1;
        end
    end
end

%% Final check, none of the dams should be flooded
if iscell(DeselectedSites)
    for i=1:numel(CheapestDam)
        if CheapestDam3(i)==1; continue; end
        DeselectedSites_final{i} = DeselectedSites{i};
    end
    
    for i=1:numel(DeselectedSites_final)
        if isrow(DeselectedSites_final{i})==0; DeselectedSites_final{i}=DeselectedSites_final{i}'; end
    end
    
    
    if exist('DeselectedSites_final')
        DeselectedSites_Final_vector = horzcat(DeselectedSites_final{:})';
        FloodCheck = sum(ismember(CheapestDam3(:),DeselectedSites_Final_vector(:)));
        fprintf('Flood check = %0.0f\n',FloodCheck);
    end
end

%% Check which index gives problem

% for i=1:numel(CheapestDam3)
%     x=sum(find(DeselectedSites_final{i}==899));
%     if x>=1
%         i
%     end;
% end

%% New COETotRD and RDPnet
COEAlls = zeros(size(COEAll));
COEAlls(COEAlls==0)=NaN;
PnetAlls = zeros(size(PnetAll));
PnetAlls(PnetAlls==0)=NaN;
ross = zeros(size(ros));
ross(ross==0)=NaN;
coss = zeros(size(cos));
coss(coss==0)=NaN;
SysIDs = zeros(size(SysID));
SysIDs(SysIDs==0)=NaN;
latss = zeros(size(lats));
latss(latss==0)=NaN;
lonss = zeros(size(lons));
lonss(lonss==0)=NaN;
LakeSurfacesAlls = zeros(size(LakeSurfacesAll));
LakeSurfacesAlls(LakeSurfacesAlls==0)=NaN;
accAlls = zeros(size(accAll));
accAlls(accAlls==0)=NaN;
BIDs = zeros(size(BIDall));
BIDs(BIDs==0)=nbasin;
CIDs = zeros(size(CIDall));
CIDs(CIDs==0)=CID;

if iscell(DeselectedSites)
    for i=1:numel(CheapestDam3)
        if CheapestDam3(i)>1
            COEAlls(CheapestDam3(i))=COEAll(CheapestDam3(i));
            PnetAlls(CheapestDam3(i))=PnetAll(CheapestDam3(i));
            
            SysIDs(CheapestDam3(i))=SysID(CheapestDam3(i));
            ross(CheapestDam3(i))=ros(CheapestDam3(i));
            coss(CheapestDam3(i))=cos(CheapestDam3(i));
            
            latss(CheapestDam3(i))=lats(CheapestDam3(i));
            lonss(CheapestDam3(i))=lons(CheapestDam3(i));
            
            LakeSurfacesAlls(CheapestDam3(i))=LakeSurfacesAll(CheapestDam3(i));
            
            accAlls(CheapestDam3(i))=accAll(CheapestDam3(i));
            BIDs(CheapestDam3(i))=BIDall(CheapestDam3(i));
            CIDs(CheapestDam3(i))=CIDall(CheapestDam3(i));
        end
    end
    
    fprintf('Total potential before deselection: %0.0f GWh\n',sum(PnetAll(~isnan(PnetAll))));
    fprintf('Total potential after deselection: %0.0f GWh\n',sum(PnetAlls(~isnan(PnetAlls))));
    
    %Deselecting dams that flood existing Grand dams
    for i=1:numel(DeselectedGrandSites_unSort)
        if isempty(DeselectedGrandSites_unSort{i})==1; continue; end
        for j=1:numel(DeselectedGrandSites_unSort{i})
            %fprintf('%d %d. %0.2f\n',i,DeselectedGrandSites_unSort{i}(j),COEAlls(i))
            COEAlls(i)=NaN; 
            PnetAlls(i)=NaN;
            SysIDs(i)=NaN;
            ross(i)=NaN;
            coss(i)=NaN;
            latss(i)=NaN;
            lonss(i)=NaN;
            LakeSurfacesAlls(i)=NaN;
            accAlls(i)=NaN;
            BIDs(i)=NaN;
            CIDs(i)=NaN;
        end
    end
    fprintf('Total potential after deselecting Grand: %0.0f GWh\n',sum(PnetAlls(~isnan(PnetAlls))));
    
    DPIDs = numel(find(SysIDs==1));
    PIDs = numel(find(SysIDs==2));
    RDIDs = numel(find(SysIDs==3));
end

%% Show dams and lake
% figure(2)
% clear img
% img{1}= truecolorsc(RDlake_Opt{CheapestDam3(1)},flipud(gray));
% for i=1:numel(CheapestDam3)
%     if i==1; continue; end
%     if CheapestDam3(i)==1;
%         img{i} = img{i-1};
%         continue;
%     end
%     
%     img{i} = burnmask(img{i-1}, ~RDlake_Opt{CheapestDam3(i)});
% end
% %
% image(img{end}); axis image;
% title('Second selection')
% hold on
% plot(co(CheapestDam3),ro(CheapestDam3),'.r','markersize',15)
% hold off
% 
% %% Show top-tot-teen graph
% % Dams can still overlap when they're in different streams
% indices_select = CheapestDam(find(~isnan(ro(CheapestDam3))));
% figure(5); clf;
% hold on;
% for i=1:numel(indices_select)
%     zbot = Z(ro(indices_select(i)),co(indices_select(i)));
%     ztop = zbot + OptDH(indices_select(i));
%     carea = acc(ro(indices_select(i)),co(indices_select(i)));
%     x = [log(carea) log(carea)];
%     y = [zbot ztop];
%     plot(x,y, 'linewidth',2,'color',rand(1,3));
% end;
% hold off;
% xlabel('log(area)');
% ylabel('Elevation');

% [lo,la] = setltln(acc, Rw, r_damst(10), c_damst(10)); %Coordinates outlets
