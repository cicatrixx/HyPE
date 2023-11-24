function val = naneq(data1, data2, roundto)
%% Check if datasets with nans are equal. Matlab does not consider nan=nan!
val=0;

if ~exist("roundto","var")
    roundto=5;
end
% round to user specified decimal points
data1=round(data1,roundto);
data2=round(data2,roundto);

% first check if all the not nan vals are equal
isnotnaneq=all(data1(~isnan(data1))==data2(~isnan(data2)),'all');

if isnotnaneq
    % check if nan vals occur in the same indexes
    if all(isnan(data1(:))==isnan(data2(:)),'all')
        val=1;
    end
end

end