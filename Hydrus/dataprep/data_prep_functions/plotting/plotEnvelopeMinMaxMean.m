function [l1,lminmax]= plotEnvelopeMinMaxMean(idata, icolor,plottype,x,swapxy,meanltype)
% idata = 2D matrix with same num of rows as the number of columns in x
% icolor = RGB color for plot
% plottype 0=meanwminmaxenvelope, 1=meanwminmaxenvelope+lines, 2=onlymean, 3=meanwminmaxenvelope+scatter
% x = 1D ROW vector for xaxis values if different from 1:nrows in idata 

% swapxy=1 uses idata as xaxis values and x as y values
% in idata, each column is a one data series
% plot takes row wise mean, max and min and plots it in given color
% https://www.mathworks.com/help/matlab/creating_plots/line-plot-with-confidence-bounds.html

if ~exist("x","var")
    x = 1:size(idata,1);
end
imean = mean(idata,2);

% create evelope
xconf = [x x(end:-1:1)] ;
yconf = [max(idata,[],2)' flipud(min(idata,[],2))'];
lminmax = [min(idata,[],2)' max(idata,[],2)'];

% swap xy data if indicated else keep xy
if exist("swapxy","var")
    xconf1=yconf;
    yconf1=xconf;
    x1=idata;
    y1=x;
    xmean=imean;
    ymean=x;
else
    xconf1=xconf;
    yconf1=yconf;
    x1=x;
    y1=idata;
    xmean=x;
    ymean=imean;
end

if ~exist("meanltype","var")
    meanltype='-';
end

% plot envelope
if ismember(plottype,[0 1 3])
    p = fill(xconf1,yconf1,'red');
    p.FaceColor = icolor; %brighten(icolor,-0.5);
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.2;
end

hold on
% plot lines or scatter
if plottype==1
    plot(x1,y1,'Color',brighten(icolor,0.3));
elseif plottype==3
    scatter(x1,y1,5,brighten(icolor,0.3));
end

% plot mean
l1=plot(xmean,ymean,'LineWidth',1.5,'LineStyle',meanltype,'Color',icolor);
