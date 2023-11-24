function addAlphaLabel(lblindex,lbltype)
% Adds label (a)','(b)','(c)'... to selected plot handle.
% lblindex: 1=(a), 26=(z)
% lbltype: place label inside or outside the plot axis box
charlbl =  compose("(%s)",('a':'z').'); % {'(a)','(b)','(c)','(d)'}
if lower(lbltype)=="inside"
    text(0.025,0.95,charlbl{lblindex},'Units','normalized','FontSize',12,'fontweight','bold')
elseif lower(lbltype)=="outside"
    text(-0.05,0.95,charlbl{lblindex},'Units','normalized','FontSize',12,'fontweight','bold')
else 
    disp("Label type assigned is not supported")
end
end
