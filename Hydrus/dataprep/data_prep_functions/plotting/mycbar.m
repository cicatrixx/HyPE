function c=mycbar(clabel,clocation)
% Create color bar w unit label in bold font and desired font location
if ~exist('clocation','var')
    clocation ='eastoutside';
end
    c=colorbar('Location',clocation);
    c.Label.String=clabel;
    c.Label.FontWeight='Bold';
    c.Label.FontSize=11;
    end