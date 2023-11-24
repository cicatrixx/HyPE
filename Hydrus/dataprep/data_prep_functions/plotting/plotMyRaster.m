function plotMyRaster(data, clabel,dname, clog)
% creates a imagesc plot where nans are greyed out with no axis grids and
% ticks. color bar has bold label and can be specified as log scale
% clog=1 log colorscale else leave blank

imagescnan(data)
c=colorbar('Location','SouthOutside');
c.Label.String=clabel;
c.Label.FontWeight='Bold';
% colormap(hot)
% freezeColors()
xticks([]); yticks([])

if exist('dname')
    title(dname,'Interpreter',  'none')
end
if exist('clog')
    set(gca, 'ColorScale', 'log')
end

%orient('landscape')
%print('G:/GitConnect/output/TheoryPot_WrapUp/TheoryPot_wSpacing.pdf','-dpdf','-fillpage')