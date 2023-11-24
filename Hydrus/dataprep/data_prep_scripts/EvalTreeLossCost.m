% Evaluate tree cover cost
% Compare the three costs
load('G:\SurfDrive\GitConnect\data\UI\data\Tree density.mat') % unit is # of trees per km2 per cell
rate=2.46; % Tree rate (USD/mound)	 
%0.06  Tree rate (USD/kg)	 
TreeCost=rate*data;
%figure;imagescnan(TreeCost)

%% Save tree cover outputs
save('G:\SurfDrive\GitConnect\data\UI\data\TreeCost.mat', 'TreeCost')
load('G:\SurfDrive\GitConnect\data\UI\data\Basin\Basin_551.mat','Rw','proj')
geotiffwrite('G:\SurfDrive\GitConnect\data\UI\data\TreeCost.tif',TreeCost,Rw,'GeoKeyDirectoryTag',proj.GeoTIFFTags.GeoKeyDirectoryTag)

%% Plot Tree cost
load('G:\SurfDrive\GitConnect\data\UI\data\Basin\Basin_551.mat','outside')
TreeCost(outside)=nan;
figure;imagescnan(TreeCost)
c=colorbar;
c.Label.String='Tree Loss Cost (USD per km2)';
c.Label.FontWeight='Bold';
