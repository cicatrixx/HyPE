function plotWaterFallDecline2Data(stack1, stack2, scenlabels, stackcolors)
%Bar plot showing how values decline over the scenarios in 2 categories of data.
%Does not work for data w increasing values
% stack1=9:-2:1;
% stack2=10:-2:2;
% stack1=flip(stack1);
% stack2=flip(stack2);
nscens=max(size(stack1));
% scenlabels=1:nscens;

diff1=abs(diff(stack1));
diff2=abs(diff(stack2));

waterfall = [stack1(1) zeros(size(diff1)) stack1(end);  %1: total1
             stack2(1) zeros(size(diff2)) stack2(end);  %2: total2
             0 stack1(2:end) 0;  %3: invisible stack1
             0 stack2(2:end) 0;  %4: invisible stack2
             0 diff1 0; % diff stack 1
             0 diff2 0]; % diff stack 2

b= bar(categorical(scenlabels,scenlabels),waterfall','stacked','FaceAlpha',0.7,'EdgeAlpha',0);

%% Make 3 and 4 stacks invisible
for i=3:4
    b(i).FaceAlpha=0;
    b(i).EdgeAlpha=0;
end

% Change color for total and succesive bars
for i=1:2
    b(i).FaceColor=stackcolors(i,:);
    b(i+4).FaceColor=brighten(stackcolors(i,:),0.2);
    b(i+4).FaceAlpha=0.3;
    
    %     b(i).EdgeColor=cl_stack(i,:);
    %b(i+2).EdgeColor=brighten(cl_stack(i,:),0.2);
end
grid on
%legend([b(1),b(2),b(5)],[strcat("Total for ",stacklabels), '\Delta change'])
legend([b(1),b(5)],{'Total';'\Delta change'})
end
