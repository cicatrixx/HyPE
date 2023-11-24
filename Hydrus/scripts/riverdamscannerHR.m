%% Riverdam scanner to determine width and height in search radius

counter=0;
for  k = 1:numel(outlets)
    % k=32;
    fprintf('%d RD dam height-width #%d of %d\n',nbasin,k,numel(outlets))
    
    if isnan(ro(k))==1; continue; end;
    
    % First load high res dem
    
    RDwinsize_dem=30;
    %RDdemHR = dem_window_corners(latOut(k), lonOut(k), RDwinsize_dem, root_bil, continent_in); %Call in function to set up high res dem window
    RDdemHR = hs3s_window('con',latOut(k), lonOut(k), RDwinsize_dem, root_bil,continent_in);
    %RDdemHR = hs3s_window('dem',latOut(k), lonOut(k), RDwinsize_dem, root_bil,continent_in);
    RDdemHR = single(RDdemHR);
    
    
    % To control for return in demHR function
    isfile_outlets(k)=0;
    if RDdemHR==0;
        isfile_outlets(k)=1;
        ro(k) = NaN;
        co(k) = NaN;
        lonOut(k)=NaN;
        latOut(k)=NaN;
        dfhdamRD{k} = 0;
        dfldamRD{k} = 0;
        ZiCenterRD(k) = 0;
        continue;
    end;
    
    [nrRD,ncRD]=size(RDdemHR);
    
    % Second find lowest point in 15s highres window
    RDwin15s=2;
    
    firstrowRD_15s = max(1,RDwinsize_dem+1-RDwin15s);
    lastrowRD_15s  = min(nrRD,RDwinsize_dem+1+RDwin15s);
    firstcolRD_15s = max(1,RDwinsize_dem+1-RDwin15s);
    lastcolRD_15s  = min(ncRD,RDwinsize_dem+1+RDwin15s);
    
    RDdemHR_15s = RDdemHR(firstrowRD_15s:lastrowRD_15s,firstcolRD_15s:lastcolRD_15s); %Select small 15s window
    [RDr_mindemHR,RDc_mindemHR]=find(RDdemHR_15s==min(RDdemHR_15s(:))); %Find lowest point
    %Transformed focal point based on demHR_15s
    %RDrh_width = RDwinsize_dem+1; %r c of focus point of large high res dem
    %RDch_width = RDwinsize_dem+1;
    RDrh_width2 = RDwinsize_dem+1+(RDr_mindemHR(1)-RDwin15s-1); % New coordinates
    RDch_width2 = RDwinsize_dem+1+(RDc_mindemHR(1)-RDwin15s-1); % In case multiple lowest point, first option
    
    ZiCenterRD(k) = RDdemHR(RDrh_width2,RDch_width2); %Corrected Zfocus
    %figure(2);clf;imagesc(RDdemHR);axis image
    % %Check
    % figure(1);clf;imagesc(RDdemHR);axis image;colormap(jet);
    % figure(1);hold on;
    % plot(RDch_width,RDrh_width,'c.','markersize',20);
    % plot(RDch_width2,RDrh_width2,'r.','markersize',20);
    % hold off
    % figure(2);clf;imagesc(RDdemHR_15s);axis image;colormap(jet);
    % figure(2);hold on; plot(RDc_mindemHR,RDr_mindemHR,'r.','markersize',20); hold off
    
    %% Height difference in dem window used as dam height
    RDdemMax(k) = max(RDdemHR(:));
    
    if MiniHydro_special==1 & Qoutlets_design(k) < bighydro_cutoff
        MiniFlag(k)=1;
        if (RDdemMax(k)-ZiCenterRD(k)-1) > MiniDamMax;
            RDZdiff(k) = MiniDamMax;
       
        else
            RDZdiff(k) = RDdemMax(k)-ZiCenterRD(k)-1;
        end
        
    else
        MiniFlag(k)=0;
        RDZdiff(k) = RDdemMax(k)-ZiCenterRD(k)-1;
    end
    
    dfhdamRD{k} = 1:RDZdiff(k)';
    
    
    %% Based on dam height we determine dam width
    
    if ASYMRD==0
        ldamRD=0;
        for wRD=1:RDZdiff(k)
            % Determine Ztop based on hdam and Zfocus
            clear RDZtop RDdx_w RDdy_w RDds_w RDhiground_w RDdhiground_w
            RDZtop = ZiCenterRD(k) + wRD;
            
            % Build higround map
            RDdx_w = repmat(abs(-RDwinsize_dem:RDwinsize_dem), [2*RDwinsize_dem+1, 1]);
            RDdy_w = RDdx_w';
            RDds_w = sqrt(RDdx_w.^2+RDdy_w.^2);
            
            % Determine minimum distance from Zfocus cell to higround
            RDhiground_w = RDdemHR>RDZtop;
            RDdhiground_w = RDds_w(RDhiground_w);
            
            ldamRD(wRD) = max(0.5,min(RDdhiground_w))*90*2;
            
            if isempty(ldamRD(wRD))==1;ldamRD(wRD)=0;end;
            
            %checkdfld
            %fprintf('lenght = %f meters.\n',ldamRD(wRD))
            %figure(1);clf;imagesc(RDhiground_w);axis image;colormap(jet);hold on
            %plot(RDch_width2,RDrh_width2,'r.','markersize',20); hold off
            
        end
        
        dfldamRD{k} = ldamRD;
        
    elseif ASYMRD==1
        
        [lrwidth{k} genwidth{k}]=AsymDamLenght(root_bil,continent_in,latOut(k),lonOut(k),RDwinsize_dem,RDZdiff(k),0);
        
        dfldamRD = lrwidth;
    end
    
end

