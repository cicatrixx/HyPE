% Evaluate maximum probable GLOF paths based on Huggel et al. 2004 Table 1
% Created By    : Sanita Dhaubanjar on 14 Nov 2020
% Created For	: SustaIndus project WP2
%=========================
clc
clear all
close all

root_data=fullfile(pwd,"data\UI\data");
load(fullfile(root_data,'UI500m_ArcGIS.mat'), 'fdir','outside','dem')
load(fullfile(root_data,'potentialGLOFs2.mat'))
glofpts=data;
cell_sz_m = 500; %in m
saveout=1;
Z=dem;

%% Loop through each GLOF risk point and identify pathway
%minslope= 0.03; %2deg from Huggel et al. 2004 Table1 Max. travel distance for lake outburst flood (flood wave) and Allen et al. 2019
%fout='G:\SurfDrive\GitConnect\data\UI\data\GLOFrisk0_03_r4';

mindeg=3; %2-3deg for flood flow, 11 deg for avalanche/debris flow from Huggel et al. 2004 Table1 Max. travel distance for lake outburst flood (flood wave) and Allen et al. 2019
minslope=round( tan(deg2rad(mindeg)),2); % 0.03-0.05 %rounding is important as paths are sensitive
fout=fullfile(root_data,sprintf('GLOFpath_%0.0fdeg_Re2',mindeg));

% get row, col of GLOF hazard lakes
[rglof, cglof]=find(glofpts>0);

figure
subplot(2,1,1)
imagescnan(maskBasin(Z,~outside))
title(sprintf('Max GLOF travel paths at min slope=%0d deg',mindeg))
hold all
plot(cglof,rglof,".k",'MarkerSize',10)

% allocate space
glofid=0*rglof;
glofpath_dist_km=0*rglof;
glofpath_vdrop_m=0*rglof;
glofpath_pathidx={};
glofpath_matrix=zeros(size(fdir));

fprintf('\n\n\n\nEvaluting max probable travel path for %d lakes at min slope=%0d deg\n',length(rglof),mindeg)
for g=1:length(rglof)
    glofid(g)=glofpts(rglof(g), cglof(g));
    [glofpath_dist_km(g), glofpath_vdrop_m(g), glofpath_pathidx{g}] = getGLOFMaxProbableTravelPath(rglof(g),cglof(g),Z,fdir,minslope,outside, cell_sz_m);
    glofpath_matrix(glofpath_pathidx{g})=glofpath_matrix(glofpath_pathidx{g})+1;
end

% Summary plot on glof paths
subplot(2,1,2)
yyaxis left
plot(glofpath_dist_km,'.-')
ylabel("Maximum hazardous path for GLOF (km)")
hold all
yyaxis right
plot(glofpath_vdrop_m,'.-')
ylabel("Head drop for GLOF (m)")
grid on

%% Plot matrix of pathways only
figure
imagescnan(maskBasin( glofpath_matrix, ~outside))
hold all
plot(cglof,rglof,".r",'MarkerSize',10)
title(sprintf('Max GLOF travel paths at min slope=%0d deg',mindeg))
text(cglof,rglof,categorical(glofid))
%% Check if glofpath_dist_km<2km
glofid(glofpath_dist_km<2)
%% Save to matlab
if saveout
    save(strcat(fout,'.mat'), 'glofid', 'glofpath_dist_km', 'glofpath_vdrop_m', 'glofpath_pathidx', 'glofpath_matrix','minslope')
    savemat2Pantpetiff(strcat(fout,'.tif'), glofpath_matrix)
    
    %% Save glof data to csv
    glofsummary=sortrows(table(glofid, glofpath_dist_km, glofpath_vdrop_m));
    writetable(glofsummary, strcat(fout,'.csv') )
    disp("Output files saved")
    
    % For saving as PDF
    % set(gcf, 'Color', 'w');
    % export_fig('G:\SurfDrive\GitConnect\data\UI\GLOFpathwaysat0_05.pdf')
    % export_fig('G:\SurfDrive\GitConnect\data\UI\GLOFpathwaysat0_03.pdf')
    
    % orient('landscape')
    % print('G:\SurfDrive\GitConnect\data\UI\CompiledDistance+CostMaps3.pdf','-dpdf','-bestfit')
end


%% Setup function getGLOFMaxProbableTravelPath
function [glofpath_dist_km, glofpath_vdrop_m, glofpath_idx] = getGLOFMaxProbableTravelPath(glof_r,glof_c,Z,fdir,minslope,outside,cell_sz_m)
%Trace downstream of glacial lake to get max probable maximum travel
%distance as per Huggel et al. 2004 Table 1 for lake outburst flood. Beyond
%these worst-case run-out distances, severe damages are not expected.
%The maximum travel distance is evaluated such that the tan(angle)>=minslope of the travel
%path from the source lake (expressed as the angle of the horizontal with
%a line from the starting point to the farthest point of deposition).


defdirs
%cell_sz_m = 500; %in m

h0=Z(glof_r,glof_c);%in m
%fprintf('!!!GLOF path search started for ID: %d!!!\n',glofid(g))

%set counter vars
slope=Inf;
dist=0; % horizontal distance
vdrop=0; %vertical drop
glofpath_idx=sub2ind(size(fdir),glof_r, glof_c ); % start at glof location itself
%% Loop through downstream cells until slope<minslope
while slope>= minslope ||   slope<0% due to DEM inconsistencies vdrop can sometimes be negative so just ignore these
    if outside(glof_r,glof_c), break; end %Stop at outlet
    
    %move downstream
    d = fdir(glof_r,glof_c);
    glof_r = glof_r+drow(d);
    glof_c = glof_c+dcol(d);
    idx=sub2ind(size(fdir),glof_r, glof_c );
    plot(glof_c,glof_r,".r")
    
    %get downstream distance based on whether it is diagonal move or not
    if d<=4 %moving E,W,N,S
        dist=dist+cell_sz_m; %in m
    else %if moving diagonaly
        dist=dist+sqrt(2*cell_sz_m^2); %in m
    end
    
    %get downstrem head and resulting slope
    vdrop=h0-Z(glof_r,glof_c);%in m
    slope=vdrop/dist;
    
    %archive
    glofpath_idx=[glofpath_idx, idx];
end
glofpath_dist_km=dist/1000; %return vals in km
glofpath_vdrop_m=vdrop;
fprintf('!!!GLOF path search stopped at %0.2fkm with head drop of %dm !!!\n',glofpath_dist_km,glofpath_vdrop_m)
end