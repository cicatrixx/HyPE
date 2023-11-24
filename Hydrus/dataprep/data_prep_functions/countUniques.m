function value_counts = countUniques(data,ntotal)
% evaluate unique values in data and list the number of occurence of each
% ignore nan cells
[C,~,ic] = unique(single(data(~isnan(data))));
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];

if exist('ntotal','var')
    value_counts(:,3)=value_counts(:,2)/ntotal*100;
end