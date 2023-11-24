 function [sumdata, ndata] = getRegionalSum(myregions, IDlist, data2sum)
% Creates groupwise sum. Similar to accumarray or groupsummary but gives
% user flexiblity to get group sums in the order they specify

% Inputs:
%    myregions - Unique group indexes or region IDs
%    IDlist -    Vector of region ID list to use for dat grouping
%    data2sum -  Vector of datavalue to sum by group. same length as IDlist
%
% Outputs:
%    output1 - Description
%    output2 - Description

for st=1:length(myregions)
    selidxs=IDlist==myregions(st);
    sumdata(st,1) = sum(data2sum(selidxs));
    ndata(st,1) = sum(selidxs);    
end
end