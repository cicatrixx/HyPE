function [basinEnergy, Hgross, idxo, idx_nbr] = evalChannelTheoreticalPot(Z_m, fdir, Q_m3s, rsi, channel, flowdist, verbose)
% Evaluate theoretical potential using channel method
% Inputs:   Z (in m), fdir (1-9), Q (m3/s), flowdist matrix of same size,
%           fdir in Patrick dirs version
%           flowdist distance from basin outlet in terms of # of cells
%           rsi river stretch interval in terms of # of cells
%           verbose is 0/1 for displaying plot and progres texts
% Outputs:  basinEnergy in GWh and Hgross in m
% Created By    : Sanita Dhaubanjar
% Created For	: SustaIndus project WP2
% Last Updated  : 30 Apr 2020
%=======================================
% Decided to keep the matrix version for consistency in the two functions
% for theoretical potential
% confirmed that traceneighbor loop goes from one outlet to another. also
% loop stays in same cell if fdir=9

disp("========evalChannelTheoreticalPot========")
% Load directional metrics for Patrick dirs
run(fullfile(pwd,'Hydrus','scripts','defdirs.m' ))

%% Channel Potential: Get sites at equal distance
% Find outlets based on hydrological distance from main outlet
% sites = channel==1 & rem(flowdist, rsi)==0; % select pixels at equal distance intervals AND within channel network
% idxo = find(sites);
% [ro,co] = ind2sub(size(Z_m), idxo);

[ro,co,idxo] = find_outletsInChannel(channel,flowdist,rsi);
fprintf("%d sites identified at cell distance = %d  \n", numel(ro),rsi)

% Plot identified sites
% figure;imagesc(Z);axis image; grid on
% hold all; plot(co,ro,'og','markersize',5) %plot outlets
% title(sprintf("%d sites identified at do of %d cells  \n", numel(ro),do))

%% Channel Potential: Get Qin at each siteHgross and trace donwstream to get Hgross
[nr,nc] = size(Z_m);
Qin = nan(size(Z_m)); %
Hgross = nan(size(Z_m)); %
idx_nbr=nan*idxo;

st=1;
for k=1:numel(ro)
    % get current cell address
    i = idxo(k);
    r = ro(k);
    c = co(k);

    % get Qin
    Qin(r,c) = Q_m3s(r,c);
    Zin = Z_m(r,c);

    % Loop through cells downstream to reach end of search distance
    for j=1:rsi
        % find downstream neighbour
        d = fdir(i);
        rj = r+drow(d);
        cj = c+dcol(d);

        % skip cell if downstream nbr is outside domain - not likely for Davids files as
        % outlets are only selected inside domain and when edge of domain
        % is reached fdir is 9 so cell is stuck, does not move down.
        if rj<1, disp("Skipped"), continue, end;
        if rj>nr, disp("Skipped"), continue, end;
        if cj<1,  disp("Skipped"), continue, end;
        if cj>nc, disp("Skipped"), continue, end;

        % if sink is reached stop looping further downstream
        if d==9
            si(st,:)=[r,c,k];
            if verbose, disp(["Sink at [r,c,k]",si(st,:)]), end
            st=st+1;
            break
        end

        % proceed downstream
        r=rj;
        c=cj;
        i=sub2ind(size(Z_m),r,c);
    end
    % Get Z and idx for downstream cell
    Zout = Z_m(r,c);
    Hgross(ro(k),co(k)) = Zin-Zout;

    %Zout1(k) = Z(r,c);
    idx_nbr(k) = i;
end

if verbose, disp("Hgross evaluated"), end
%% Channel Potential: Evaluate potential
g           = 9.8;          %Gravitational acceleration (m/s2)
rho         = 1000;         %Density of water (kg/m3)
hr_yr      = 8760;          %Hrs in year... (24*365 = 8760 hr/yr)
eta        = 1;             %Efficiency turbine

% Grid based results
basinPot = rho * g * eta .* Hgross.* Qin;  % power in W
basinEnergy = basinPot* hr_yr *1e-9; % energy in GWh

fprintf("Total +ve Energy GWh in cells in basin: %0.2f\n", sum(basinEnergy(basinEnergy>=0)))
fprintf("Total -ve Energy GWh in cells in basin: %0.2f\n", sum(basinEnergy(basinEnergy<0)))

%% Column based results -- Diff if I make this int
% Zin1 = single(Z(idxo)); %
% Qin1 = single(Q(idxo)); %
% Zout1 = single(Z(idx_nbr)); %
% basinPot2 = rho * g * eta .* (Zin1-Zout1).* Qin1 ;  % power in W
% basinEnergy2 = basinPot2* hr_yr ; % energy in Wh
% disp("Channel pot evaluated")

%% Outputs: Plot distribution of channel potential
if verbose
    nbasinchannel=numel(ro);
    nvalid=sum(Hgross(:)>=0);
    fprintf("Num sites in basin: %d\n",nbasinchannel)
    fprintf("Num -ve sites in basin: %d (%0.2f%% of channel cells)\n", sum(Hgross(:)<0),sum(Hgross(:)<0)/nbasinchannel*100)
    fprintf("Total +ve Hgross in cells in basin: %0.2f\n", sum(Hgross(Hgross>=0)))
    fprintf("Total -ve Hgross in cells in basin: %0.2f\n", sum(Hgross(Hgross<0)))

    mybasin=Z_m>0; %
    % Get not -ve values from matrix
    basinEnergy_channel = basinEnergy;
    basinEnergy_channel(basinEnergy<0)=0;

    figure;
    subplot(1,2,1)
    colormap(gray)
    imagesc(mybasin);axis image;grid on % plot gray background
    hold on;
    [rch,cch]=ind2sub(size(Z_m), find(channel));
    scatter(co,ro,10+round(basinEnergy_channel(idxo)*1e-2),'MarkerEdgeColor','#D95319') %,'filled') % the 1e-1 multiplilcation is only for scaling symbol size - 1 necessary as 0 is not allowed
    plot(cch,rch,'.b','markersize',1) %plot channel
    colordata = colormap;
    colordata(1,:) = .8*[1 1 1];
    colormap(colordata);
    legend(sprintf("HP: %0.3g- %0.0f GWh",min(basinEnergy_channel(basinEnergy_channel>0)),max(basinEnergy_channel(:))), "Channel")
    title(sprintf("Channel Potenial for rsi = %d (Total: %0.2f GWh at %d valid sites)", rsi,nansum(basinEnergy_channel(:)),nvalid))
    %h=colorbar('Location','eastoutside');
    h.Label.String='MW';

    subplot(1,2,2)
    bar(sort(basinEnergy_channel(~isnan(basinEnergy_channel)))) % plot only not nans
    grid on
    set(gca, 'YScale', 'log')
    ylabel("Theoretical Potential at Channel Cell (GWh)")
    title('Distribution of Potential Energy in the Channel ')
end
disp("========evalChannelTheoreticalPot========")
end