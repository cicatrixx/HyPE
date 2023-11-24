%% Create the Rw Georeference structure

% Order of ID based on aggregation.m Continents = {'AFR';'ASIA';'AUS';'CAM';'EUR';'NAM';'SAM'};
if strcmp('NAM',continent_in)
    R = georasterref();
    R.RasterSize = [nr nc];
    R.Latlim = [25 60];
    R.Lonlim= [-138 -52];
    R.ColumnsStartFrom = 'north';
    CID=6;
    % EUR
elseif strcmp('EUR',continent_in)
    R = georasterref(); % Geo ref raster
    R.RasterSize = [nr nc];
    R.Latlim = [12 62];
    R.Lonlim= [-14 70];
    R.ColumnsStartFrom = 'north';
    CID=5;
    % SAM
elseif strcmp('SAM',continent_in)
    R = georasterref(); % Geo ref raster
    R.RasterSize = [nr nc];
    R.Latlim = [-56 15];
    R.Lonlim= [-93 -32];
    R.ColumnsStartFrom = 'north';
    CID=7;
    % AFR
elseif strcmp('AFR',continent_in)
    R = georasterref(); % Geo ref raster
    R.RasterSize = [nr nc];
    R.Latlim = [-35 38];
    R.Lonlim= [-19 55];
    R.ColumnsStartFrom = 'north';
    CID=1;
    % ASIA
elseif strcmp('ASIA',continent_in)
    R = georasterref(); % Geo ref raster
    R.RasterSize = [nr nc];
    R.Latlim = [-12 61];
    R.Lonlim= [57 180];
    R.ColumnsStartFrom = 'north';
    CID=2;
    % AUS
elseif strcmp('AUS',continent_in)
    R = georasterref(); % Geo ref raster
    R.RasterSize = [nr nc];
    R.Latlim = [-56 -10];
    R.Lonlim= [112 180];
    R.ColumnsStartFrom = 'north';
    CID=3;
    % CAM
elseif strcmp('CAM',continent_in)
    R = georasterref(); % Geo ref raster
    R.RasterSize = [nr nc];
    R.Latlim = [5 39];
    R.Lonlim= [-119 -60];
    R.ColumnsStartFrom = 'north';
    CID=4;
end

[lu_lat, lu_lon] = setltln(basin, R, firstrow, firstcol); % Handy function to calc lat lon from row
[rb_lat, rb_lon] = setltln(basin, R, lastrow, lastcol);

Rw = georasterref(); % Geo ref raster for basin
Rw.RasterSize = [nrw ncw];
Rw.Latlim = [rb_lat lu_lat];
Rw.Lonlim= [lu_lon rb_lon];
Rw.ColumnsStartFrom = 'north';

%% Alternative method withou mapping_tool box
%     %NAM
% if strcmp('NAM',continent_in)
%     Rlat = linspace(60,25,nr);
%     Rlon = linspace(-138,-52,nc);
%     CID=6;
%     % EUR
% elseif strcmp('EUR',continent_in)
%     Rlat = linspace(62,12,nr);
%     Rlon = linspace(-14,-70,nc);
%     CID=5;
%     % SAM
% elseif strcmp('SAM',continent_in)
%     Rlat = linspace(15,-56,nr);
%     Rlon = linspace(-93,-32,nc);
%     CID=7;
%     % AFR
% elseif strcmp('AFR',continent_in)
%     Rlat = linspace(38,-35,nr);
%     Rlon = linspace(-19,55,nc);
%     CID=1;
%     % ASIA
% elseif strcmp('ASIA',continent_in)
%     Rlat = linspace(61,-12,nr);
%     Rlon = linspace(57,180,nc);
%     CID=2;
%     % AUS
% elseif strcmp('AUS',continent_in)
%     Rlat = linspace(-10,-56,nr);
%     Rlon = linspace(112,180,nc);
%     CID=3;
%     % CAM
% elseif strcmp('CAM',continent_in)
%     Rlat = linspace(39,5,nr);
%     Rlon = linspace(-119,-60,nc);
%     CID=4;
% end
% 
% Rlatw = linspace(Rlat(firstrow),Rlat(lastrow),nrw);
% Rlonw = linspace(Rlon(firstcol),Rlon(lastcol),ncw);

% latOutn = [   46.2629
%    48.1007
%    48.2090
%    48.1882
%    48.2090
%    48.0090
%    47.7756
%    47.7173
%    47.4173
%    47.3381
%    45.9879
%    45.8712
%    51.4177
%    51.3802
%    51.9053
%    51.8845
%    51.4969
%    51.8094
%    45.4795
%    45.4545
%    45.4545
%    45.9754
%    45.5962
%    46.6422
%    46.2379
%    46.3046
%    46.5046
%    49.6258
%    49.6425
%    49.6509
%    49.6383
%    49.6258
%    47.8798
%    47.1381
%    47.0089
%    47.8215
%    47.8548
%    47.8715
%    42.5125
%    43.4418
%    43.4126
%    42.6208
%    43.3959
%    43.2084];
% 
% lonOutn = [ -117.1344
%  -118.9720
%  -119.0803
%  -119.0595
%  -119.0803
%  -118.8803
%  -118.6470
%  -118.5886
%  -118.2886
%  -118.2095
%  -116.8594
%  -116.7427
%  -122.2888
%  -122.2513
%  -122.7763
%  -122.7555
%  -122.3680
%  -122.6805
%  -116.3510
%  -116.3260
%  -116.3260
%  -116.8469
%  -116.4677
%  -117.5136
%  -117.1094
%  -117.1761
%  -117.3761
%  -120.4971
%  -120.5137
%  -120.5221
%  -120.5096
%  -120.4971
%  -118.7512
%  -118.0094
%  -117.8803
%  -118.6928
%  -118.7261
%  -118.7428
%  -113.3842
%  -114.3134
%  -114.2843
%  -113.4926
%  -114.2676
%  -114.0801];
% 
% %% Krijg ik met deze lat lons dezelfde row cols als de andere methode?
% 
% for j=1:numel(ro)
%     if isnan(ro(j))==1; ron(j)=NaN; con(j)=NaN; continue; end;
%     
%     [ron(j),con(j)] = setpostn(Q,Rw,latOutn(j),lonOut(j));
% end
% 
% isequal(ro,ron)
% isequal(co,con)
% 
% %%
% figure(1);clf;imagesc(log(Q));axis image;colormap(flipud(gray));hold on
% plot(co,ro,'.r','markersize',20)
% plot(con,ron,'.b','markersize',15)
% hold off
%
% %% ANTWOORD: Net niet helemaal
