function mycmap=mycbrewer(cname,nmap)
% Many cmaps in Cbrewer start w white as the first color which is useless
% for TS data. this function cushions the cmap to remove the lightest
% shades

tmp=cbrewer2(cname,nmap+2);
mycmap=tmp(3:end,:);