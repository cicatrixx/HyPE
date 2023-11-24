%% Skip if max flowdist is < do (distance between outlets)

if max(flowdist(:)) < do;
    disp('Flowdist is less than distance between outlets');
    save(matfile,'-v7.3','do');   
    returnID=1;
    return
end

%% Find outlets and create windows

[ro,co] = find_outlets_dis(flowdist,acc,do);
outIdx = sub2ind(size(acc),ro,co);

%Make variables
for k=1:numel(ro)
    Accoutlets(k)=0;
    PopDisplacedOpt(k)=0;
    Quakerate(k)=0;
    COETotRD(k)=0;
    RDPnet(k)=0;
    Qoutlets_design(k)=0;
    Qoutlets_design_LF(k)=0;
    RDRegion_id(k)=0;
    RDCountry_id(k)=0;
    latOut(k)=0;
    lonOut(k)=0;
    OptDH(k)=0;
    OptDL(k)=0;
    OptPop(k)=0;
    OptLV(k)=0;
    RDDepth(k)=0;
    OptInv(k)=0;
    RDP(k)=0;
    RDlakeSurface{k}=[];
    RDVolumeLake{k}=[];
    RDVolumeLake15s(k)=0;
    RDSurfaceLake15s(k)=0;
    Zoutlets(k)=0;
    DisOutlet(k)=0;
    OptSpecCap(k)=0;
    PopDisplacedOpt(k)=0;
    Zoutlets(k)=0;
    Qoutlets(k)=0;
    Qoutlets_design(k)=0;
    Qoutlets_design_LF(k)=0;
    RDRegion_id(k)=0;
    RDCountry_id(k)=0;
    RDDepth(k)=0;
    Accoutlets(k)=0;
    
end

% figure(1);clf;imagesc(log(Q));colormap(flipud(gray));axis image; hold on
% plot(co,ro,'r.','markersize',15); hold off
% sum(isnan(ro))
% figure(2);ax2=subplot(1,3,2);imagesc(log(Q));colormap(flipud(gray));axis image; hold on
% plot(co,ro,'r.','markersize',15); hold off
%% Q-maps with ot without water consumption

if waterconsumption==1
    Q=Qwc;
    Qdesign=Qdesign_wc;
    Qdesign_mean=Qdesign_mean_wc;
    Qdesign_LF=Qdesign_LF_wc;
end

%% Lat lon coordinates outlets

for  k = 1:numel(ro)
    [latOut(k), lonOut(k)] = setltln(acc, Rw, ro(k), co(k)); %Coordinates outlets
    %     latOut(k) = Rlatw(ro(k))-((Rlatw(2)-Rlatw(1))/2);
    %     lonOut(k) = Rlonw(ro(k))-((Rlonw(2)-Rlonw(1))/2);
    Depth_Rivermouth(k) = D(ro(k),co(k)); % Depth of locations from outlet
    flowdist_rivermouth(k) = flowdist(ro(k),co(k)); % Distance from outlet
    Qoutlets(k) = Q(ro(k),co(k));
    outlets(k) = sub2ind(size(Z),ro(k),co(k));
    Accoutlets(k) = acc(ro(k),co(k));
end

%% Deselect outlet and - No dam in river mouth

if rivermouth_constr==1
    if numel(ro)<outletdeselect;
        disp('To few ro and co');
        save(matfile,'-v7.3','do');  
        returnID=1;
        return
    end
end

%% Rivermouth constraint
if rivermouth_constr==1
    
    [~,sortQidx]=sort(Accoutlets,'descend');
    
    %First two for sure
    ro(sortQidx(1:outletdeselect))=NaN;
    co(sortQidx(1:outletdeselect))=NaN;
%     ro(sortQidx(50:(end)))=NaN;
%     co(sortQidx(50:(end)))=NaN;
    
    %This method causes additional deselection after lakes in the river
    %system, especially in Congo and Amazon
    %Second deselect dams on mainstream xkm inland and on rivers with a certain depth
    for k=1:numel(ro)
        if flowdist_rivermouth(k) < rivermouth_inland & Depth_Rivermouth(k) > depth_cutoff
        ro(k)=NaN;
        co(k)=NaN;
        end
    end
    
    %Other 2 maybe depending on river depth
%     if strcmp('SAM',continent_in)
%         if nbasin==1
%             if numel(sortQidx) < 8; RDi = numel(sortQidx); else RDi=33; end
%             for i=3:RDi
%                 if Depth_Rivermouth(sortQidx(i)) > depth_cutoff
%                     ro(sortQidx(i))=NaN;
%                     co(sortQidx(i))=NaN;
%                 end
%             end
%         end
%     end
end
% sum(isnan(ro))

% figure(1); hold on
% plot(co,ro,'b.','markersize',15); hold off
% figure(2); subplot(1,3,2); hold on
% plot(co,ro,'b.','markersize',15); hold off
% fd=flowdist;fd(Q>100)=0;cmap=jet(256);cmap(1,:)=[1 1 1];
% figure(2);ax3=subplot(1,3,3); imagesc(fd);colormap(cmap);axis image;

%% Deselect locations in WDPA regions

if protarea_constr==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        if WDPA_PL10(ro(k),co(k))== 10 || WDPA_PL20(ro(k),co(k))==20
            ro(k)=NaN;
            co(k)=NaN;
        else
        end
    end
else
end

%% Navigation constraint for navigability (4m)

if navi_constr==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        if D(ro(k),co(k))> depth_cutoff
            ro(k)=NaN;
            co(k)=NaN;
        else
        end
    end
else
end
% [rd, cd] = setpostn(Q, Rw, -3.128930, -60.027996);
% figure(2); hold on
% plot(cd,rd,'g.','markersize',25); 
% plot(co,ro,'b.','markersize',15); hold off
% sum(isnan(ro))
%% Mangrove constraint

if mangrove_constr==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        if MC(ro(k),co(k))==1
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
end

%% Mainstream-constraint - No dams on basin mainstream

if mainstream_constr==1
    [~,MainIdx] = max(Q(:));
    mainstrm = find_mainstream(Q,MainIdx,adir);
    outlets(ismember(outlets,mainstrm))= NaN;
    ro(find(isnan(outlets)))=NaN;
    co(find(isnan(outlets)))=NaN;
end

%% seismicnoaa-constraint - No dams in earthquake areas

if seismicnoaa_constr==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        if seismicnoaa(ro(k),co(k))==1
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
end

%% MiniHydro deselection

if MiniHydro_deselect==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        if Q(ro(k),co(k))< bighydro_cutoff
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
end

%% MiniHydro Selection

if MiniHydro_select==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        if Q(ro(k),co(k))> bighydro_cutoff
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
end

%% Existing dams de-selection

if ExistingDams_constraint==1
    for k= 1:numel(ro)
        if isnan(ro(k))==1; continue; end;
        if NoDamsLand(ro(k),co(k))== 1
            ro(k)=NaN;
            co(k)=NaN;
        end
    end
end

%% No dams before first GrandDam on mainstream

if BeforeFirstDam_constraint==1
    if r_damst(1)~=99;
        NoBeforeFirstMainstrDam = FirstMainstrmDam(Q,adir,acc,r_damst,c_damst,GrandIdx);
        outlets(ismember(outlets,NoBeforeFirstMainstrDam))= NaN;
        ro(find(isnan(outlets)))=NaN;
        co(find(isnan(outlets)))=NaN;
    end
end

%% Deselect dam location with zero discharge value (can happen due to water consumption)

for k= 1:numel(ro)
    if isnan(ro(k))==1; continue; end;
    if Q(ro(k),co(k))== 0
        ro(k)=NaN;
        co(k)=NaN;
    end
end


%% if all NaN break loop

if sum(isnan(ro))==numel(ro);
    disp('First all NaN break');
    save(matfile,'-v7.3','do');
    returnID=1;
    return
end

%% Transform
for k=1:numel(ro)
    
    if isnan(ro(k))==1; continue; end;
    
    Zoutlets(k) = Z(ro(k),co(k));
    Qoutlets(k) = Q(ro(k),co(k));
    Qoutlets_design(k) = Qdesign(ro(k),co(k));
    Qoutlets_design_LF(k) = Qdesign_LF(ro(k),co(k));
    RDRegion_id(k) = Regions(ro(k),co(k));
    RDCountry_id(k) = Countries(ro(k),co(k));
    RDDepth(k) = D(ro(k),co(k));
    DisOutlet(k) = Dis(ro(k),co(k)); % Distance to powerline (km)
end

%% Slackflow for ecological reasons
if slackflow_constraint==1
    Qoutlets_design = Qoutlets_design *( 1-slackflow);
end

%% Control for NaNs in Dis map (slightly inconsistent with other maps because different sea map was used)
for k=1:numel(DisOutlet)
    if isnan(DisOutlet(k))==1; ro(k)=NaN; co(k)=NaN; end;
end

% figure(2); hold on
% plot(co,ro,'b.','markersize',15); hold off

%% Seismic PIK
% 5 percent cost increase at seismic hazard > 4

if seismicpikcost==1
    for k=1:numel(ro)
        if isnan(ro(k))==1; Quakerate(k)=NaN; continue; end;
        
        if seismicpik(ro(k),co(k)) > 4
            Quakerate(k) = 0.05;
        end
    end
else
    for k=1:numel(ro)
        Quakerate(k) = 0;
    end
end

%% Report deselection
no=numel(outlets);
noNAN=sum(isnan(ro));
fprintf('Number of outlets: %d.\n', no);
fprintf('Number of deselected: %d.\n', noNAN);

%% Create windows for Qdecreaser

for k = 1:numel(outlets)
    
    if isnan(ro(k))==1; continue; end;
    
    inlet_win = dowin;
    ZiWinfr(k) = ro(k)-inlet_win;
    ZiWinlr(k) = ro(k)+inlet_win;
    ZiWinfc(k) = co(k)-inlet_win;
    ZiWinlc(k) = co(k)+inlet_win;
    
    if ZiWinfr(k)<1 || ZiWinlr(k)>nrw || ZiWinfc(k)<1 || ZiWinlc(k)>ncw
        ro(k)=NaN;
        co(k)=NaN;
    else
    end
    
    if isnan(ro(k))==1; continue; end;
    
    
    fdir_inlet_win{k}      = fdir(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    adir_inlet_win{k}      = adir(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    Q_inlet_win{k}         = Q(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    acc_inlet_win{k}       = acc(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    Z_inlet_win{k}         = Z(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    Pop_inlet_win{k}       = Pop(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    WDPA_PL10_inlet_win{k} = WDPA_PL10(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    WDPA_PL20_inlet_win{k} = WDPA_PL20(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    LandValue_inlet_win{k} = LandValue(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    Dis_inlet_win{k}       = Dis(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    flowdist_inlet_win{k}  = flowdist(ZiWinfr(k):ZiWinlr(k),ZiWinfc(k):ZiWinlc(k));
    
    
end

WinProb=sum(isnan(ro))- noNAN;
fprintf('Number of window problems: %d.\n', WinProb);

index_inlet_win = sub2ind([(2*inlet_win+1) (2*inlet_win+1)],inlet_win+1,inlet_win+1);

%% if all NaN break loop

if sum(isnan(ro))==numel(ro);
    disp('all NaN brak loop');
    save(matfile,'-v7.3','do');
    save(matfileCOEPOT,'-v7.3','do');    
    returnID=1;
    return
end

% figure(2); hold on
% plot(co,ro,'b.','markersize',15); hold off

%%
% disp('Saving Depth')
% matpath = fullfile(root, sprintf('output\\%s\\%d', continent_in, nbasin));
% matfile = fullfile(root, sprintf('output\\%s\\%d\\Basin%d_Depth.mat', continent_in, nbasin, nbasin));
% if ~isdir(matpath)
%     mkdir(matpath);
% end
%
% save(matfile,'-v7.3','RDDepth');
%
% continue

% disp('Test');
% save(matfile,'-v7.3','do','ro','co','ro_org','co_org');
% returnID=1;
% return