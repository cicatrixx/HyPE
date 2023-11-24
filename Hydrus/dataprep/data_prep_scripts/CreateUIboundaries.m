%% Create catchment boundaries for UI at 5km and 500m and merge them. Generate 500m accuflux and ldds in same extent as 5km.

%%% At 5km, new UI boundaries based on Sonu's SPHY model outputs, skipping
% downstream subbasin. there were some trailing cells so to manually generated
% mask in arcgis to clean these up!
%%% At 500m, snapped Sonu's UI outlets to river based on 500m accuflux and
% ldd in arcgis. Used new outlets to generate 500m catchments. But found that
% endoheric basin in north east is parsed as part of Indus main so manually
% changed ldd in Qgis Serval. Used altered ldd to generate altered accuflux
% which is used to generate final 500m catchments. Found that extent of ldd
% file was incorrect so reran to clip extents for ldd and outlets file.
%%% Merged 5km and 500m catchments by first creating unique IDs and summing
% them. 99 and 9900 used for subcatchments outside the UI
% Convert all map files to be used in later stages to tif

% USES GDAL_TRANSLATE!
% USES PC_RASTER!

%%
close all
clear all
root5km="G:/SPHY/Sonu/full_runs_sanita/data";
root500m="G:/SPHY/Arthur/500m";
%cd(root500m);

process5km=10;
process500m=10;
correctldd500m=0;
resample5kmto500m=0;
merge5km500m=10;

%% Create 5km mask
if process5km
    my5km_outlets=1;  %Not used anymore as it generates 0s where cat 6,7,9 and 11 were. Better to remove these from Sonu's catchments instead
    if my5km_outlets
        %% Create new outlet points only in UI
        % cd(root);
        ifile=fullfile(root5km,"outlets5.map");
        ofile=fullfile(root5km,"outlets5_UIonly.map");
        convertmap2tif(ofile); % Final outlet5km
        
        %system(strcat("aguila ",ifile))
        command=strcat("pcrcalc ",ofile,"= if(",ifile," != 6 and ",ifile," != 7 and ", ifile," != 9 and ",ifile," != 11, ",ifile,")");
        system(command)
        
        system(strcat("aguila ",ofile))
        %% Generate new UI catchments from new outlet points
        lddfile=fullfile(root5km,"ldd5.map");
        newldd=ofile;
        catfile=fullfile(root5km,"UIcatchment5k.map");
        system(strcat("pcrcalc ",catfile," = catchment(",lddfile,", ",newldd,")"))
        
        system(strcat("aguila ",catfile))
    end
    
    %% Generate new UI catchments from Sonu's catchments
    ifile=fullfile(root5km,"catchments5.map");
    ofile=fullfile(root5km,"catchments5_UIonly.map");
    
    system(strcat("aguila ",ifile))
    command=strcat("pcrcalc ",ofile,"= if(",ifile," != 6 and ",ifile," != 7 and ", ifile," != 9 and ",ifile," != 11, ",ifile,")");
    system(command)
    
    system(strcat("aguila ",ofile))
    
    %% Translate manual mask to clean up sonu's basin from  .geotiff to .map
    ifile=fullfile(root5km,'catchments5_UI_clmask.tif');  % Mask created in Arcmap after manually cleaning up the 0 grid cells in the LI
    mask = strrep(ifile,".tif","matlabGDAL.map");
    command = sprintf("gdal_translate -of PCRaster %s %s", ifile,mask);
    
    system(command)
    
    %% Clean up basin mask from Sonu's catchments
    ifile=fullfile(root5km,"catchments5_UIonly.map");
    ofile=fullfile(root5km,"UIcatchment_5km_cl.map");
    
    command=strcat("pcrcalc ",ofile,"= if(!",mask,", ",ifile,")");
    system(command)
    
    % Change basin ID 0 to 99 for easier tracking
    system(strcat("pcrcalc ",ofile," = if(",ofile,"==0, 99,",ofile,")"))
    system(strcat("aguila ",ofile))
    
    convertmap2tif(ofile); % Final subbasinmask5km
    
    
    %% Create mask
    finalmask=fullfile(root5km,"mask5km.map");
    system(strcat("pcrcalc ",finalmask," = boolean(",ofile,")"))
    system(strcat("aguila ",finalmask))
    
    convertmap2tif(finalmask); % Final basinmask5km
    disp("Prepared 5km outlets, subbasins and basin mask")
    
end

%% Convert 5km .map to .tif
convertmap2tif(fullfile(root5km,'accuflux.map'));    %Final acc5km
convertmap2tif(fullfile(root5km,'dem.map'));         %Final dem5km
convertmap2tif(fullfile(root5km,'ldd5.map'));         %Final ldd5km
%% Create 500m mask

if process500m
    %% Correct 500m ldd extents and accuflux for endoheric basin
    if correctldd500m
        %% Convert 500m ldd raster from QGIS .tif to .asc -- manually clip to correct extents
        lddtif=fullfile(root500m,"ldd_500m_v2SD.tif");  % Corrected in Qgis manually
        newlddtif = fullfile(root500m,"ldd_500m_v3SD.tif") ;
        tmpasc = strrep(lddtif,".tif","tmp.asc");
        
        %command = sprintf("gdal_translate %s %s", newlddtif,tmpasc);
        system(sprintf("gdal_translate -projwin -1670000.0 840000.0 -230000.0 30000.0 %s %s", lddtif, newlddtif)) %Final ldd500km
        system(sprintf("gdal_translate %s %s", newlddtif, tmpasc))
        
        %% Use asc2map for ldd conversion
        newldd = strrep(newlddtif,".tif",".map");
        extent=fullfile(root500m,"extent500m.map");  % from arthur
        command = sprintf("asc2map --clone %s -a -L %s %s", extent, tmpasc, newldd); %convert to LDD type
        system(command)
        %system(strcat("aguila ", newldd))
        
        % Repair ldd
        system(sprintf("pcrcalc %s = lddrepair(%s)", newldd, newldd))
        
        %% create new accuflux
        newacc=fullfile(root500m,"accuflux_500m_v3SD.map");
        system(strcat("pcrcalc ",newacc,"= accuflux(", newldd,", 1)"));
        system(strcat("aguila ", newacc))
        
        %% convert new accuflux to tiff
        convertmap2tif(newacc);  %Final accuflux500km
        disp("Prepared corrected 500m accuflux")
        
    end
    
    %% Convert 500m outlet raster from ArcGIS .tif to .asc to .map becoz it doesnt work otherwise
    outlettif2map=1;
    if outlettif2map
        ifile=fullfile(root500m,"outlets5_500m.tif");  % Mask created in Arcmap manually
        ifiletrim=fullfile(root500m,"outlets5_500m_v3SD.tif");  % Mask created in Arcmap manually
        system(sprintf("gdal_translate -projwin -1670000.0 840000.0 -230000.0 30000.0 %s %s", ifile, ifiletrim))  %Final outlet500m
        
        tmpfile = strrep(ifiletrim,".tif","tmp.asc");
        command = sprintf("gdal_translate %s %s", ifiletrim,tmpfile);
        system(command)
        
        ofile = strrep(ifiletrim,".tif",".map");
        command = sprintf("gdal_translate %s %s", tmpfile,ofile);
        system(command)
        %system(strcat("aguila ", ofile))
    end
    
    %% Generate 500m UI catchments
    lddfile=fullfile(root500m,"ldd_500m_v3SD.map");
    ptfile=fullfile(root500m,"outlets5_500m_v3SD.map");
    catfile=fullfile(root500m,"UIcatchment_500m_v3SD.map");
    system(strcat("pcrcalc ", catfile," = catchment(",lddfile,", ",ptfile,")"))
    
    system(strcat("aguila ", catfile))
    
    % Clean up 500m catchments
    cat500m=strrep(catfile,"_v3","_cl_v3");
    
    % Change basin ID 0 to 99 for easier tracking
    system(strcat("pcrcalc ",cat500m," = if(",catfile,"==0, 99,",catfile,")"))
    system(strcat("aguila ",cat500m))
    convertmap2tif(cat500m); % Final subbasinmask500m
    
    
    %% Create 500m mask
    mask500m=fullfile(root500m,"mask500m_v3SD.map");
    system(sprintf("pcrcalc %s = boolean(%s)",mask500m,cat500m))
    system(strcat("aguila ",mask500m))
    disp("Prepared 500m catchments")
    convertmap2tif(mask500m); % Final basinmask500m
    
end

%% Resample 5km catchment to resolution and extent of 500m catchment - SLOWEST STEP
cat500m=fullfile(root500m,"UIcatchment_500m_cl_v3SD.map");
newcat5kmTo500m=fullfile(root500m,"UIcatchment_5kmTo500m_cl2_v3SD.map");
if resample5kmto500m
    cat5km=fullfile(root5km,"UIcatchment_5km_cl.map");
    cat5kmTo500m=fullfile(root500m,"UIcatchment_5kmTo500m_cl_v3SD.map");
    system(strcat("resample --clone ", cat500m, " ", cat5km, " ", cat5kmTo500m))
    system(strcat("aguila ",cat5kmTo500m))
    
    % Multiply 5km basin ID by 1000
    system(sprintf("pcrcalc %s = scalar(%s)*100", newcat5kmTo500m, cat5kmTo500m))
    
    % Create mask 5kmto500m
    mask5kmTo500m=fullfile(root500m,"mask5kmTo500m_v3SD.map");
    system(sprintf("pcrcalc %s = boolean(%s)", mask5kmTo500m, cat5kmTo500m))
    system(strcat("aguila ",mask5kmTo500m))
    
end

%% Merge 5km and 500m catchments
if merge5km500m
    % Get intersection of not nan values from 500m and 5km
    intersectmap=fullfile(root500m,"UIcatchment_intersect5km500m_v3SD.map");
    unionmap=fullfile(root500m,"UIcatchment_union5km500m_v3SD.map");
    system(sprintf("pcrcalc %s = scalar(%s) + %s",intersectmap, cat500m, newcat5kmTo500m)) % This generates intersection only
    system(strcat("aguila ",intersectmap))
    
    % Add nan values to intersection:
    % Create raster w non-nan values from inputs. For cells w non-nan vals in
    % multiple input, take the first. so first take vals from intersection of
    % 500m and 5km then take remainders from 500m
    system(sprintf("pcrcalc %s = nominal(cover(%s, %s, scalar(%s)))",unionmap, intersectmap, newcat5kmTo500m, cat500m))
    system(strcat("aguila ",unionmap))
    
    % Set 99 from 500m to nan
    unionmap2=fullfile(root500m,"UIcatchment_union5km500m_cl_v3SD.map");
    system(sprintf("pcrcalc %s = if(%s!=99,%s)",unionmap2,unionmap,unionmap))
    system(strcat("aguila ",unionmap2))
    
    % Create union mask
    unionmask5km_500m=fullfile(root500m,"mask_union5km500m_v3SD.map");
    system(sprintf("pcrcalc %s = boolean(%s)", unionmask5km_500m, unionmap2))
    system(strcat("aguila ",unionmask5km_500m))
    
    %% Translate mergefiles from .map to .geotiff
    convertmap2tif(unionmap2); % Final unioned subbasin500m
    convertmap2tif(unionmask5km_500m); % Final unioned basinmask500m
    
    disp("Merged 5km and 500m catchments")
    
end

%% Fix extents of dem and convert to tif
olddem=fullfile(root500m,"dem_500m.map");
demfixed=fullfile(root500m,"dem_500m_v3SD.map");
extent=fullfile(root500m,"extent500m.map");  % from arthur
system(strcat("resample --clone ", extent, " ", olddem, " ", demfixed))
system(strcat("aguila ",demfixed))
convertmap2tif(demfixed); % Final dem500m

%% Test scripts gives the same only not nan vals of the two
% system("pcrcalc tmp3.map = mask500m.map or mask5kmTo500m.map")
% system("pcrcalc tmp3.map = mask500m.map and mask5kmTo500m.map")
%  system("pcrcalc UIcatchment_intersect5km500m.map = boolean(UIcatchment_500m_cl.map) and boolean(UIcatchment_5kmTo500m_cl.map)")
% system("pcrcalc tmpUnion2.map = nominal(if(mask5kmTo500m.map and mask500m.map, scalar(UIcatchment_500m_cl.map) + UIcatchment_5kmTo500m_cl2.map, if(mask5kmTo500m.map, UIcatchment_5kmTo500m_cl2.map, if(mask500m.map, scalar(UIcatchment_500m_cl.map)))))")
% system("pcrcalc tmpUnion3.map = if(mask500m.map, UIcatchment_500m_cl.map)")
% system("pcrcalc tmpUnion4.map = if(mask500m.map, if(mask5kmTo500m.map then scalar(UIcatchment_500m_cl.map) + UIcatchment_5kmTo500m_cl2.map))")
% system("pcrcalc tmpUnion5.map = if(mask500m.map then if(mask5kmTo500m.map then scalar(UIcatchment_500m_cl.map) + UIcatchment_5kmTo500m_cl2.map else scalar(UIcatchment_500m_cl.map)) else 0)")
% system("pcrcalc tmpUnion6.map = if(mask5kmTo500m.map then UIcatchment_500m_cl2.map + UIcatchment_5kmTo500m_cl2.map else UIcatchment_500m_cl2.map)")
% get union

%% Fix extents of accuflux
% newaccfixed=fullfile(root500m_BB,"accuflux_500m_v3SD.map");
% extent=fullfile(root500m_BB,"extent500m.map");  % from arthur
% system(strcat("resample --clone ", extent, " ", newacc, " ", newaccfixed))
% system(strcat("aguila ",newaccfixed))

%% Test scripts to find best way to clip ldd
% root500m_BB="G:/SPHY/Arthur/500m/getBasinBounds";
%
% system("resample --clone extent500m.map outlets5_500m.map outlets5_500m_v2SD.map")
%
% % i) clip first in QGIS and then asc2map and lddrepair
% system(sprintf("cd %s",root500m_BB))
% system("gdal_translate -projwin -1670000.0 840000.0 -230000.0 30000.0 ldd_500m_v2SD.tif ldd_500m_v2SD_clipped.asc")
% system("asc2map --clone extent500m.map -a -L ldd_500m_v2SD_clipped.asc ldd_500m_v2SD_clipped.map")
% system("pcrcalc ldd_500m_v2SD_clipped2.map = lddrepair(ldd_500m_v2SD_clipped.map)")
% system("pcrcalc diffldd2.map=scalar(ldd_500m_v2SD_clipped2.map)-scalar(ldd_500m_v2SD_clipped.map)")
% %system("aguila diffldd2.map")
% system("pcrcalc cat_clipped2.map=catchment(ldd_500m_v2SD_clipped2.map, outlets5_500m_v2SD.map)")
% system("aguila cat_clipped2.map")
%
% % ii) asc2map to nominal and then resample and lddrepair
% system("gdal_translate ldd_500m_v2SD.tif ldd_500m_v2SD.asc")
% system("asc2map --clone ldd_500m.map -a -N ldd_500m_v2SD.asc ldd_500m_v2SD.map")
% system("resample --clone extent500m.map ldd_500m_v2SD.map ldd_500m_v2SD_clippedT2.map")
% system("pcrcalc ldd_500m_v2SD_clippedT22.map = lddrepair(ldd(ldd_500m_v2SD_clippedT2.map))")
% system("pcrcalc difflddT22.map=scalar(ldd_500m_v2SD_clippedT22.map)-scalar(ldd_500m_v2SD_clippedT2.map)")
% %system("aguila diffldd2.map")
% system("pcrcalc cat_clippedT22.map=catchment(ldd_500m_v2SD_clippedT22.map, outlets5_500m_v2SD.map)")
% system("aguila cat_clippedT22.map")
%
% % compare the two methods
% system("pcrcalc difflddT1_T2.map=scalar(ldd_500m_v2SD_clippedT22.map)-scalar(ldd_500m_v2SD_clipped2.map)")
% system("pcrcalc cat_T1_T2.map=scalar(cat_clippedT22.map)-scalar(cat_clipped2.map)")
%
%compare w arthur's ldd
system("pcrcalc difflddT1_Arthur.map=scalar(ldd_500m_v2SD_clipped2.map)-scalar(ldd_500m_correct_extent_sound.map)")

