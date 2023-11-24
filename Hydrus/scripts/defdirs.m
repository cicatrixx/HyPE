% defdirs.m : define directions
% Small script that ensures consistent use of constants throughout scripts
% and functions

% version 1.0,  12 oct 2006,  P.W.Bogaart
% version 1.1,  24 jan 2013

%% Define directional constants
if exist('NORTH','var'), return, end;   % Only run when needed
%disp('Defining directions');

% Define 4 cardinal directions
NORTH = 1;
WEST  = 2;
EAST  = 3;
SOUTH = 4;

cardinal_dirs = 1:4;

% Define 4 diagonal directions.
% These are optional, that's why they are defined later
NWEST = 5;
NEAST = 6;
SWEST = 7;
SEAST = 8;

SELF  = 9;

diagonal_dirs = 5:8;

% Define a special value to denote sink points (depressions in the DEM)
% where flow direcetion in undefined
SINK        = 99;

% Define directional constants as usded by ESRI
ESRIdirs = [64, 16, 1, 4, 32, 128, 8, 2];

% Define opposite directions
antidir = zeros(1,9);
antidir(NORTH) = SOUTH;
antidir(SOUTH) = NORTH;
antidir(EAST)  = WEST;
antidir(WEST)  = EAST;

antidir(NWEST) = SEAST;
antidir(NEAST) = SWEST;
antidir(SWEST) = NEAST;
antidir(SEAST) = NWEST;

antidir(SELF) = SELF;

% Define the relative inter-node/cell distances
reldelta = zeros(1,9);
reldelta(1:4) = 1.0;
reldelta(5:8) = sqrt(2.0);
reldelta(9) = 0.0;

% Define shift in row/columns
drow = zeros(1,9);
drow(NWEST) = -1;
drow(NORTH) = -1;
drow(NEAST) = -1;
drow(WEST)  =  0;
drow(SELF)  =  0;
drow(EAST)  =  0;
drow(SWEST) =  1;
drow(SOUTH) =  1;
drow(SEAST) =  1;

dcol = zeros(1,9);
dcol(NWEST) = -1;
dcol(NORTH) =  0;
dcol(NEAST) =  1;
dcol(WEST)  = -1;
dcol(SELF)  =  0;
dcol(EAST)  =  1;
dcol(SWEST) = -1;
dcol(SOUTH) =  0;
dcol(SEAST) =  1;

%eof

%% Define hydropower constants
minQ=0.1; %minimum annual avg flow constrain
