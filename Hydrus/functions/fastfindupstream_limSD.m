function basin = fastfindupstream_limSD(acc,fdir,drow,dcol, outlet,showplot,cntchk)
%% fast-find catchment upstream of given outlet position
% Changes:
    % counterchk as input
    % added evalcell
    % added showplot
    

if ~exist('cntchk','var')
    cntchk = 1000;
end

[nr,nc] = size(acc);

[~,order] = sort(acc(:),'descend'); %plot(acc(order)) 

basin = zeros(nr,nc, 'logical'); % reset at empty map
evalcell=zeros(nr,nc, 'logical'); % tracker to check if cell is evaluated

basin(outlet) = true; % outlet is lowest grid cell in catchment
N = nr*nc;
firstk = find(order==outlet); % position within 'order' where outlet is.
pold=0;
counter=0;

for k = firstk : N
    % use counter to limit number of cells upstream of outlet chosen
    counter=counter+1;
    if counter==cntchk
       if showplot; fprintf("Counter 1000 at order: %d \n",order(k)), end
       continue
    end
    
    i = order(k);
    % dont do anything at outlet location
    if basin(i), continue, end
    
    r = rem(i-1, nr)+1;
    c = fix((i-1)/nr)+1;
    % find neighbour using fdir
    d = fdir(i);
    rj = r+drow(d);
    cj = c+dcol(d);
    if rj<1,  continue, end
    if rj>nr, continue, end
    if cj<1,  continue, end
    if cj>nc, continue, end
    
    % copy catchment membership status from downstream nbr
    basin(r,c) = basin(rj,cj);
    evalcell(r,c)=1;
    
    % report
    p = fix(100*k/N);
    if showplot && p>pold
        fprintf('Progress: %d%%.\n', p);
        pold=p;
    end
end
if showplot
    imagescnan(evalcell);axis image
    title(sprintf("# of cells evaluated: %0.0f",sum(evalcell(:))))
    fprintf("Counter ended at order: %d \n",order(k))
end