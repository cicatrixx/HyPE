function [basinEnergy, Hgross] = evalCell2CellTheoreticalPot(Z_m, fdir, Q_m3s, verbose)
% Evaluate theoretical potential using cell by cell method
% Inputs:   Z (in m), fdir and Q (m3/s) matrix of same size,
%           fdir in Patrick dirs version
%           verbose is 0/1 for displaying plot and progres texts
%           !!! nan cells in Z is checked to ID only basin cells to loop
%           through !!!
% Outputs:  basinEnergy in GWh and Hgross in m
% Created By    : Sanita Dhaubanjar
% Created For	: SustaIndus project WP2
% Last Updated  : 30 Apr 2020
%=======================================

disp("========evalCell2CellTheoreticalPot========")
%% Theoretical Potential: Set up variables and params
g           = 9.8;          %Gravitational acceleration (m/s2)
rho         = 1000;         %Density of water (kg/m3)
hr_yr      = 8760;          %Hrs in year... (24*365 = 8760 hr/yr)
eta        = 1;             %Efficiency turbine

% Load directional metrics moving according to Patrick dirs
run('defdirs.m' )

% Initialize placeholders
[nr,nc] = size(Z_m);
Hgross = nan(size(Z_m)); %
Hout = nan(size(Z_m)); %

% Counter
pold=0; %counter

% Get index for not nan cells to read through
order = find(~isnan(Z_m));
N = numel(order);

%% Theoretical Potential: Get Hgross at each cell
% Loop through all cells, ID downstream cell using fdir and get elevation out
for k = 1 : N
    i = order(k);
    
    % get outlet rc
    r = rem(i-1, nr)+1;
    c = fix((i-1)/nr)+1;
    
    % find downstream neighbour
    d = fdir(i);
    rj = r+drow(d);
    cj = c+dcol(d);
    
    % skip cell if
    if rj<1,  continue, end;
    if rj>nr, continue, end;
    if cj<1,  continue, end;
    if cj>nc, continue, end;
    
    % get headloss at inlet
    Zin = Z_m(r,c);
    Zout= Z_m(rj,cj);
    if Zin<Zout
        %keyboard
    end
    Hout(r,c) = Z_m(rj,cj);
    Hgross(r,c) = Zin-Zout;%For Gross head (m)
    
    % report
    p = fix(100*k/N);
    if verbose && p>pold
        fprintf('Progress: %d%%.\n', p);
        pold=p;
    end;
end
if verbose, disp("Hgross evaluated"), end

%% Theoretical Potential: Get basin potential
basinPower = rho * g * eta .* Hgross .* Q_m3s ;  % power in W
basinEnergy = basinPower* hr_yr *1e-9; % energy in GWh
fprintf("Total +ve Energy GWh in cells in basin: %0.2f\n", sum(basinEnergy(basinEnergy>=0)))
fprintf("Total -ve Energy GWh in cells in basin: %0.2f\n", sum(basinEnergy(basinEnergy<0)))

if verbose
    nbasincell=sum(~isnan(Hgross(:)));
    nvalid=sum(Hgross(:)>=0);
    fprintf("Num cells in basin: %d\n",nbasincell)
    fprintf("Num -ve cells in basin: %d (%0.2f%% of total)\n", sum(Hgross(:)<0),sum(Hgross(:)<0)/nbasincell*100)
    fprintf("Total +ve Hgross in cells in basin: %0.2f\n", sum(Hgross(Hgross>=0)))
    fprintf("Total -ve Hgross in cells in basin: %0.2f\n", sum(Hgross(Hgross<0)))
    
    basinEnergy_cell=basinEnergy;
    basinEnergy_cell(basinEnergy<0)=0;%%
    figure;
    subplot(1,2,1)
    imagescnan(basinEnergy_cell);axis image;grid on
    title(sprintf("Cell by Cell Potenial (Total: %0.2f GWh at %0.3g%% valid basin cells)",nansum(basinEnergy_cell(:)), nvalid/nbasincell*100))
    h=colorbar('Location','eastoutside');
    h.Label.String='GWh';
    set(gca,'ColorScale','log')
    % Color nans as grey
%     colordata = colormap;
%     colordata(1,:) = .95*[1 1 1];
%     colormap(colordata);
%     freezeColors
    %
    tmp=sort(basinEnergy_cell(:));
    subplot(1,2,2)
    bar(tmp(tmp>0))
    grid on
    set(gca, 'YScale', 'log')
    ylabel("Theoretical Potential at Cell (GWh)")
    title('Distribution of Potential Energy in Cell ')
end
disp("========evalCell2CellTheoreticalPot========")
end