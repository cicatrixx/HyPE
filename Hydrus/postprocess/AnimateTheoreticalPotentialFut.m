
%Create final plot for: dQ vs dHP and Z vs HP
% other plots are created in SubBasinWisePotential_Fut script
clear all
close all
addpath(genpath('Hydrus'))
run('myVarNames_Fut.m')
saveastif=0;
minQ=0.1; %in m3/s
miscCheck=0;
plotfigs=0;
% Load theoretical outputs
load(fullfile(rootof,"FutRuns_Figs_Analysis","fig_theory","FutTheoreticalPot_totals.mat"),...
    'fname','dorange','channel_GWh','Qannual_m3s_archive','channel_basinEnergy') % Q and HP data are both already re-ordered and masked

% Load subbasin boundaries
datapath=fullfile(pwd,"data","UI","data");
load(fullfile(datapath,'UI500m_ArcGIS.mat'),'catchments','basinlabels','outside','dem','channel')

disp("Loaded theory potential files")


%% Read correct order of models
modorder=readtable(fullfile(rootof,"FutRuns_Figs_Analysis","FutScen_Tracker.xlsx"),Sheet="GCMdetails", Range="A1:f13");
suffix=strcat("_LTavgs_",strrep(tframe,'-','_'));
myforder=[strcat(modorder.matlabname,suffix(1)); strcat(modorder.matlabname,suffix(2))];
newforder=zeros(1,length(myforder));

for f=1:length(myforder)
    newforder(f)=find(fname==myforder(f));
end
% keep historical at start
newforder=[1 newforder];
newfnames=fname(newforder');
channel_GWh_n=channel_GWh(:,newforder);
channel_basinEnergy_n=channel_basinEnergy(:,:,newforder);
Qannual_m3s_archive_n=Qannual_m3s_archive(:,:,newforder);
% prepare cleaned up names
ssp_gcm_names=extractBefore(newfnames(2:end),"_LTavgs_");
ssp_gcm_tf_names=["Historical"; strcat(ssp_gcm_names(1:12), "_",tframe(1)); strcat(ssp_gcm_names(13:24), "_",tframe(2))];
disp("Reorederd future scenarios")

%% Save theorypot layers as tif
if saveastif
    fpath2Rw = fullfile(rootf,sprintf('\\data\\%s\\Basin_UIB\\PantpeBasin_%d.mat',continent_in,nbasin));
    prefix="";
    for scen =1:size(data3D_in,3)
        if scen==2
            prefix="Mid_";
        elseif scen==14
            prefix="Far_";
        end
        tiffoutpath = fullfile(rootoffut,"TheoryPotTiffs", strcat(prefix, strrep(strrep(histrcpcornernames{scen},'.','_'),':', '_'),".tif"));
        savemat2Pantpetiff(tiffoutpath, data3D_in(:,:,scen), fpath2Rw)
    end
end


%% Setup for spatial maps
data3D_in = channel_basinEnergy_n;

% Eval prct change
data3D_in_prctchange=(data3D_in-data3D_in(:,:,1))./data3D_in(:,:,1)*100;
animate3DMatrix2GIF(data3D_in_prctchange,'Scenario - ',tfhistrcpcornernames,'tmp.GIF')

minscal=min(data3D_in,[],'all');
maxscal=max(data3D_in,[],'all');

%[rch,cch]=find(channel);

minmarker=3;
maxmarker=1800;
cdem=colormap("gray");
mlaborange=[0.8500 0.3250 0.0980];

plotprefix="Scenario";
ptitle=tfhistrcpcornernames;

%% Plot maps
for scen=1:5
    % Get scen data
    sel_basinEnergy=data3D_in_prctchange(:,:,scen); %sel the do=1 case
    idxo=  find(~isnan(sel_basinEnergy));
    potval=sel_basinEnergy(idxo);
    [ro,co]= ind2sub(size(sel_basinEnergy),idxo  );
    nvalid=numel(ro);


    %% FINAL: Plot DEM + vizualized DB + theoretical potential
    figure(scen)
    imagescnan(maskBasin(dem,~outside))
    mycbar('Elevation (m)')
    colormap(cdem(150:end,:));
    hold all

    scatter(co,ro,rescale(potval,minmarker,maxmarker,'InputMax',maxscal,'InputMin',minscal),mlaborange,'filled','MarkerFaceAlpha',0.3,'DisplayName','Theoretical Potential') %,'MarkerFaceColor','#D95319')cmap(1,:)
    %legend('Color',cdem(length(cdem)-10,:))
    title([plotprefix, ptitle{scen}], Interpreter="none")

end