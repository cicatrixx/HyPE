function regionID=getRegionID(in_r, in_c, regionMatrix)
%getRegionID - Returns region ID for not nan row, column indices provided
% r,c pair is used to generated single index for the regionMatrix
%
% Syntax:  regionID=getRegionID(in_r, in_c, regionMatrix)
%
% Inputs:
%    in_r - Vector of row indices
%    in_c - Vector of column indices. Same length as in_r
%    regionMatrix - 2D matrix with region IDs
%
% Outputs:
%    regionID - Vector with region IDs for (in_r, in_c) pairs. Same length as in_r.

%------------- BEGIN CODE --------------
regionID=nan(size(in_r));

% Get index only for not nan (R, C)
valid_idxs=find(~isnan(in_r));
rsel = in_r(valid_idxs);
csel = in_c(valid_idxs);
isel=sub2ind(size(regionMatrix), rsel, csel);

% Get region ID from region matrix for selected indexes
regionID(valid_idxs)=regionMatrix(isel);

%------------- END OF CODE --------------
