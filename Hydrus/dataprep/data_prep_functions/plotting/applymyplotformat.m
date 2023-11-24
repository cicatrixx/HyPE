function applymyplotformat(mytitle,cmap)

title(mytitle,'Interpreter','none')
if exist('cmap')
    colororder(cmap)
end
    
box on
grid on % apply major grid
%grid minor
%mygraylines=.5*[1 1 1];  % for axes borders and boxes
%set(gca,'XColor',mygraylines,'yColor',mygraylines)
set(gca,'FontName',    'Segoi UI','FontSize',11);
end
