function [rida cida ida]= MaxIDBlockFinder(acc,scalingf)
% First find highest acc idx in coarser blocks defined by scalingf, then
% return only unique ida and its r,c

%scalingf = 120; %scaling factor 15s map to 0.5deg as .5*60*60s/15s
scalingf = 10; %scaling factor 15s map to 0.5deg as .5*60*60s/15s

%% Select 0.5x0.5 degree block
[nr nc] = size(acc);

rb = linspace(scalingf,nr,nr/scalingf);
cb = linspace(scalingf,nc,nc/scalingf);

ct=0;
r1=1;
c1=1;
i=0;
clear accBlock
for r=rb
    i=i+1;
    if i==1; r1(i)=1; else r1(i)=rb(i-1)+1; end
    r2(i)=r;
    j=0;
    
    for c=cb
        j=j+1;
        if j==1; c1(j)=1; else c1(j)=cb(j-1)+1; end
        c2(j)=c;
        
        ct=ct+1;
        accBlock{ct} = acc(r1(i):r2(i),c1(j):c2(j));
        if isempty(find(accBlock{ct}))==1; id(ct)=1; rid(ct)=1; cid(ct)=1;continue;end
        [~, id(ct)] = max(accBlock{ct}(:));                         %index of highest cell in acc block
        [rid(ct), cid(ct)] = ind2sub(size(accBlock{ct}),id(ct));    %row col of index
        
    end
end


%% check
for x=1:numel(rid)
    %         x=175;
    %         figure(1);clf;
    %         ax1=subplot(1,2,1);
    %         imagesc(log(accBlock{x}));axis image;colormap(flipud(gray));hold on
    %         plot(cid(x),rid(x),'.r','markersize',15)
    %         hold off
    
    nrs=floor(x/numel(cb));
    ncs=x-(nrs*numel(cb));
    
    if nrs==0; accr1=1; else accr1=rb(nrs)+1; end; %r1 of block
    accr2=accr1+scalingf-1;                        %r2 of block
    
    if ncs==1|ncs==0; accc1=1; else accc1=cb(ncs-1)+1; end; %c1 of block
    accc2=accc1+scalingf-1;                          %c2 of block
    
    %         fprintf('%d\n',nrs);
    %         fprintf('%d\n',ncs);
    %         fprintf('%d\n',accr1);
    %         fprintf('%d\n',accr2);
    %         fprintf('%d\n',accc1);
    %         fprintf('%d\n\n',accc2);
    %
    %         accBlockCheck = acc(accr1:accr2,accc1:accc2); %acc block
    %         ax2=subplot(1,2,2);
    %         imagesc(log(accBlockCheck));axis image;colormap(flipud(gray));hold on
    %         plot(cid(x),rid(x),'.r','markersize',20)
    %         hold off
    %         linkaxes([ax1 ax2])
    
    rida(x) = rid(x)+accr1-1; %row number of big raster
    if rida(x)>nr; rida(x)=nr; end;
    cida(x) = cid(x)+accc1-1; %col number of big raster
    
    idx(x) = sub2ind(size(acc),rida(x),cida(x));
    
    %     figure(3);clf;imagesc(log(acc));axis image;colormap(flipud(gray));hold on
    %     plot(cida,rida,'.r','markersize',20)
    %     hold off
    %
    %         fprintf('r/c AccBlock = %0.2f\n',accBlock{x}(rid(x),cid(x)))
    %         fprintf('r/c Acc      = %0.2f\n\n',acc(rida(x),cida(x)))
end

%% Select only unique Idas
[idxUn, ia, ~] = unique(idx);
rida = rida(ia);
cida = cida(ia);
ida = idxUn;


rida=rida((~isnan(ida)));
cida=cida((~isnan(ida)));
ida=ida((~isnan(ida)));
end