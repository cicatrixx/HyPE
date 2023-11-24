function basin = fastfindupstream_Dis(acc,fdir,Dis,drow,dcol,outlet,dowin)
%% fast-find catchment upstream of given outlet position
[nr,nc] = size(acc);
[~,order] = sort(acc(:),'descend');

Dis_order = Dis(order);

basin = logical(zeros(nr,nc, 'single')); % reset at empty map
basin(outlet) = true; % outlet is lowest grid cell in catchment
N = nr*nc;
firstk = find(order==outlet); % position within 'order' where outlet is.
StartDis = Dis_order(firstk); % Start position
pold=0;
for k = firstk : N
    i = order(k);
    if Dis(i)-StartDis > dowin, continue, end;
    %if k==firstk, assert(i==outlet); end;
    if basin(i), continue, end; % happens at outlet location
    r = rem(i-1, nr)+1;
    c = fix((i-1)/nr)+1;
    % find neighbour
    d = fdir(i);
    rj = r+drow(d);
    cj = c+dcol(d);
    if rj<1,  continue, end;
    if rj>nr, continue, end;
    if cj<1,  continue, end;
    if cj>nc, continue, end;
    % copy catchment membership status from downstream nbr
    basin(r,c) = basin(rj,cj);
    

    
    % report
    p = fix(100*k/N);
    if p>pold
        %fprintf('fastfindupstream_Dis Progress: %d%%.\n', p);
        pold=p;
    end;
end
