%% RiverDam lake to calculate population displacement

for k=1:numel(outlets)
    fprintf('%d RD Lake Pop-Landvalue #%d of %d\n',nbasin,k,numel(outlets))
    
    if isnan(ro(k))==1; continue; end;

    %Calculate upstream cells
    RDupstream = fastfindupstream_lim(acc,fdir,drow,dcol,outIdx(k));
    
    % define box around upstream cells to speed up
    ZupstreamRD = RDupstream==1;
    %ZupstreamRD = RDupstream==1 & Z > Z(ro(k),co(k)) & Z < (Z(ro(k),co(k))+dfhdamRD{k}(end));
    cs = find(sum(ZupstreamRD)>0); % columns
    rs = find(sum(ZupstreamRD')>0); % rows
    colwin=0;
    
    fcs = max(1,cs(1) - colwin);
    lcs = min(ncw,cs(end) + colwin);
    frs = max(1,rs(1) - colwin);
    lrs = min(nrw,rs(end) + colwin);
    Zupstream_win = Z(frs:lrs,fcs:lcs);
    RDupstream_win = RDupstream(frs:lrs,fcs:lcs);
    Popupstream_win = Pop(frs:lrs,fcs:lcs);
    LandValueupstream_win = LandValue(frs:lrs,fcs:lcs);
    
    %Determine which cells are flooded
    PopRDlake=0;
    MoneyRDlake=0;
    LandValueRDlakee=0;
    clear RDlake2 hd
    
    for j=1:numel(dfhdamRD{k})
        hd=single(dfhdamRD{k}(j));
        
        RDlake2 = RDupstream_win & Zupstream_win < (Zoutlets(k)-RDDepth(k)+OptDH(k)); % Based on 15s DEM map, minus river depth + dam height based on 3s DEM
        
        %Surface
        RDlakeSurface{k}(j) = sum(RDlake2(:)==1)* (450*(450*cosd(latOut(k)))); % m2 
        
        %Volume
        Z_upstream = Zupstream_win(RDupstream_win);
        dz_h = Zoutlets(k)-RDDepth(k)+hd-Z_upstream;
        dz_h = max(0, dz_h);
        RDVolumeLake{k}(j)  = sum(dz_h * 450 * (450*cosd(latOut(k)))); % m3
        
        %Population
        RDlakeidx = find(RDlake2);
        PopRDlake(j) = sum(Popupstream_win(RDlakeidx));
        
        %Land value check with runsettings
        if landval_set==1
            LandValueRDlakee(j) = sum(LandValueupstream_win(RDlakeidx))*0.45^2 * 9; %LandValue is $/km2. Area of 15s cell is 450m^2 x 9.09 (infinite discounted yearly revenue with discount factor 10%, see perpetuity_check.m)
        else
            LandValueRDlakee(j) = 0;
        end
        
        %% Find GDP values
        %GDPSwiss = 55000; % $/cap/year
        %GDPUSA = 54000; % $/cap/year
        [GDPr,~] = find(ISOGDP(:,1)==RDCountry_id(k));
        GDPpc = ISOGDP(GDPr,3);
        if isempty(GDPr)==1; GDPpc = mean(ISOGDP(:,3)); end %if ISO value is not found, world average
        
        MoneyRDlake(j) = 2 * GDPpc * PopRDlake(j); %
    end
    
    PopDisplaced{k}=PopRDlake;
    PopCost{k}=MoneyRDlake;
    %RDLake_map{k}=RDlake2;
    LandValueRDlake{k}=LandValueRDlakee;
    
end

%% Check lake

% figure(1);clf;
% img1= truecolorsc(Zupstream_win,flipud(gray));
% % img2 = burnmask(img1, ~RDLake_map{14}{40});
% img2 = burnmask(img1, ~RDlake2);
% %             img2 = burnmask(img1, ~RDLake_map{k}{200});
% image(img2); axis image;
% title('Discharge','FontSize',14);
% %     axis([75 125 75 125]);
% set(gca,'xtick',[],'ytick',[])
% hold on;
% % plot(co(23),ro(23),'r.','markersize',15);
% hold off

%%
% [nr nc] = size(Zupstream_win)
% for r=1:nr
%     for c=1:nc
%         if RDupstream_win(r,c)==1
%             Zup(r,c) = Zupstream_win(r,c);
%         else
%             Zup(r,c) = NaN;
%         end
%     end
% end
% figure(3);clf;imagesc(Zup); axis image
% 
% %%
% for r=1:nr
%     for c=1:nc
%         Zupm(r,c) = 710 - Zup(r,c);
%     end
% end
% figure(4);clf;imagesc(Zupm); axis image
% 
% %%
% Zupmm = max(0, Zupm);
% figure(5);clf;imagesc(Zupmm); axis image
