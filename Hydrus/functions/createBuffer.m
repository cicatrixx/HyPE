function bufferedmat=createBuffer(inmat,buffersize)
% Create a new mat with buffer around cells w value 1 in inmat
% inmat is a 0/1 matrix
% buffersize in terms of no. of cells

[nr, nc]=size(inmat);
% find non zero cells
[cr, cc]=ind2sub(size(inmat),find(inmat));

bufferedmat=zeros(size(inmat),'logical');
% Loop through all non zero cells in inmat
for ri= 1:length(cr)
    r=cr(ri);
    c=cc(ri);
    % Generate buffer indices such that values >=1 or <=nc/nr, within domain
    bufferedmat(max(1,r-buffersize):min(nr, r+buffersize),max(1,c-buffersize):min(nc,c+buffersize))=1;
end
end
%bufferonly=bufferedmat-inmat;
%newchannel=inmat+createBuffer(inmat,1)+createBuffer(inmat,2);
