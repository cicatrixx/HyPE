function inlets = inlet_finderSD(Z,Qyravg,acc,fdir,WDPA_PL10,WDPA_PL20,protarea_constr,nr,nc,Zout,ZPin,deselect_mainstream,index_inlet_win,adir,flowdist,sradius)
% Return matrix the same size as input data window with 1=inlet or 0
% Inlet identified as upstream cells within search radius that have Q>=minQm3/s,
% at input elevation level and optionally outside protected areas 

% Zout=%index of center dam outlet location
% ZPin= min Z threshold for current inlet search

% Z=Z_inlet_win{k};
% Q=Q_inlet_win{k};
% acc=acc_inlet_win{k};
% fdir=fdir_inlet_win{k};
% WDPA_PL10=WDPA_PL10_inlet_win{k};
% WDPA_PL20=WDPA_PL20_inlet_win{k};
% protarea_constr=1;
% nr=nr_ZiW;
% nc=nc_ZiW;
% Zout;
% ZPin=Zin(l);

defdirs

% %Pinletsw = inlet_finder(Z=Z_inlet_win{k},
% Q=Q_inlet_win{k},
% acc=acc_inlet_win{k},
% fdir=fdir_inlet_win{k},
% WDPA_PL10=WDPA_PL10_inlet_win{k},
% WDPA_PL20 = WDPA_PL20_inlet_win{k},
% protarea_constr = protarea_constr,
% nr_ZiW,
% nc_ZiW,
% Zout,
% ZPin(l),
% deselect_mainstream,
% index_inlet_win,
% adir_inlet_win{k},
% flowdist_inlet_win{k},
% sradius);

%upstream = fastfindupstream(acc,fdir,drow,dcol,Zout);
%upstream = fastfindupstream_Dis(acc,fdir,flowdist,drow,dcol,Zout,sradius); %sradius search radius 50 cells = 25km
%upstream = fastfindupstream_DisSD5(acc,fdir,adir,flowdist,drow,dcol,Zout,0,sradius); %no plot %dowin=sradius #fastfindupstream_DisSD5(acc,fdir,adir,flowdist,drow,dcol,outlet,showplot, dowin)
upstream = fastfindupstream_DisSD4(acc,fdir,adir,flowdist,drow,dcol,Zout,0,sradius); %no plot %dowin=sradius #fastfindupstream_DisSD5(acc,fdir,adir,flowdist,drow,dcol,outlet,showplot, dowin)
 

inlets = logical(zeros(nr,nc, 'single'));

% Selection of inlets based on settings
if protarea_constr==1
    locs = find(upstream & Z>=ZPin & Qyravg>=minQ & WDPA_PL10==0 & WDPA_PL20==0)';
else
    locs = find(upstream & Z>=ZPin & Qyravg>=minQ)';
end

if deselect_mainstream==1
    mainstream = find_mainstream(Qyravg,index_inlet_win,adir);
    locs(ismember(locs,mainstream))= NaN;
    locs = locs(find(~isnan(locs)));
end

for i = locs
    r = rem(i-1, nr)+1;
    c = fix((i-1)/nr)+1;
    % find downstream nbr
    rj = r+drow(fdir(i));
    cj = c+dcol(fdir(i));
    inlets(r,c) = Z(rj,cj) < ZPin;
end

end