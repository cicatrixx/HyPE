function [x, y, idx]=bresenham(x1, y1, x2, y2, matsize)
% Code by Aaron Wetzler (2021). Bresenham optimized for Matlab
% (https://www.mathworks.com/matlabcentral/fileexchange/28190-bresenham-optimized-for-matlab),
% MATLAB Central File Exchange. Retrieved August 25, 2021.
% This code does not use any for loops and takes advantage of Matlabs
% internally optimized routines to produce a fast, optimized version of
% Bresenham's line drawing algorithm

%Format:
%               [x y]=bham(x1,y1,x2,y2)
%
%Input:
%               (x1,y1): Start position
%               (x2,y2): End position
%
%Output:
%               x y: the line coordinates from (x1,y1) to (x2,y2)
%
%Usage example:
%               [x y]=bham(1,1, 10,-5);
%               plot(x,y,'or');

x1=round(x1); x2=round(x2);
y1=round(y1); y2=round(y2);
dx=abs(x2-x1);
dy=abs(y2-y1);
steep=abs(dy)>abs(dx);

if steep
    t=dx;dx=dy;dy=t;
end

%The main algorithm goes here.
if dy==0
    q=zeros(dx+1,1);
else
    q=[0; diff(mod([floor(dx/2):-dy:-dy*dx+floor(dx/2)]',dx))>=0]; % select increment amount for next point
end
%and ends here.

if steep
    if y1<=y2
        y=[y1:y2]';
    else
        y=[y1:-1:y2]';
    end
    if x1<=x2
        x=x1+cumsum(q);
    else
        x=x1-cumsum(q);
    end
else
    if x1<=x2
        x=[x1:x2]';
    else
        x=[x1:-1:x2]';
    end
    if y1<=y2
        y=y1+cumsum(q);
    else
        y=y1-cumsum(q);
    end
end
if exist('matsize', 'var')
 idx = sub2ind(matsize, y, x);
end