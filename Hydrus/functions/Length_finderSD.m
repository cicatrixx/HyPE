function L = Length_finderSD(acc,outlet_idx,inlet_idx,latOut,Zout,Zin,coordsys, cellsz)
% Evaluate distance between inlet and outlet indexes using pythagoras theorem.
% Matrix acc is used to convert indexes to row/column values. Latitude correction
% is applied for coordsys="geographic". Cellsz is used to convert cell distance to
% nominal distances in same units as cellsz.

%Row and colum of outlet
[rOutlet, cOutlet] =  ind2sub(size(acc),outlet_idx);

%Row and colum of inlets
[rInlet, cInlet] =  ind2sub(size(acc),inlet_idx);

%First horizontal distance
if strcmp(coordsys, 'geographic')
    %get lat lon for the inlet
    a = ((rOutlet - rInlet) * cellsz) * cos(deg2rad(latOut)); % Longitudial length, corrected for latitude
elseif strcmp(coordsys, 'planar')
    a = ((rOutlet - rInlet) * cellsz); % Longitudial length
end
b = (cOutlet - cInlet)  * cellsz; % Latitudial length
hL = sqrt(a^2+b^2); %Horizontal distance

%Second correcting for vertical distance
L = sqrt(hL^2 + (Zout - Zin)^2); %Total pipelength

end