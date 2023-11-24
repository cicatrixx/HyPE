load('G:\SurfDrive\HPmodel\data\ASIA\Basin_UIB\PantpeBasin_103.mat', 'Heritages_natural', 'outside')
nat_km2=sum(Heritages_natural,'all')*.5*.5;
UIB_km2=sum(~outside,'all')*.5*.5;
fprintf("Area of natural heritages: %0.2f sq km",nat_km2)
fprintf("Area of UIB: %0.2f sq km",UIB_km2)
fprintf("Percent of UIB area covered by nat heritages: %0.2f %%",nat_km2/UIB_km2*100)