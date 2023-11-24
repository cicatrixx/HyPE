function addHPcategoryLines(addtextlabel,VorH,xpos)
% Add Vertical (VorH=="v") or horizontal (VorH=="h") lines separating HP size wise categories

% Techlim=0.5;% $/kWh
% Econlim=0.1;% $/kWh

% % from Hoes et al
% HPclass=["Mega(>1000 MW)", "Large (10-1000 MW)", "Small (1-10 MW)", "Mini (0.1-1 MW)", "Micro (0.005-0.1 MW)", "Pico (<0.005 MW)"];
% HPsz_greaterthan=[1000, 10, 1, 0.1, 0.005]/1000*365*24; %MW converted to GWh

% add lines and labels for different sizes of HP -- from Siddiqui
HPclass=["Mega (>1000 MW)", "Large (500-1000 MW)", "Medium (50-500 MW)", "Small (5-50 MW)", "Mini (0.15-5 MW)", "Micro (0.005-0.15 MW)", "Pico (<0.005 MW)"];
HPsz_greaterthan=[1000	500	50	10	5	0.005]/1000*365*24; %MW converted to GWh assuming year round production

if VorH=="h"
    %add lines for different categories of HP
    for i=1:length(HPclass)-1
        yline(HPsz_greaterthan(i),'--k');
    end

    if addtextlabel
        if ~exist('xpos','var')
            xmax=xlim();
            xpos=xmax(2)*.95;
        end

        text(repmat(xpos,length(HPclass),1),[HPsz_greaterthan*1.5 HPsz_greaterthan(end)*.5] ,HPclass,'FontAngle','italic','DisplayName','')%'FontWeight','bold',
    end
elseif VorH=="v"
    %add lines for different categories of HP
    for i=1:length(HPclass)-1
        xline(HPsz_greaterthan(i),'--k');
    end

    if addtextlabel
        ymax=ylim();
        text([HPsz_greaterthan HPsz_greaterthan(end)*.5],repmat(ymax(2)*1.1,length(HPclass),1) ,HPclass,'FontWeight','bold','FontAngle','italic','DisplayName','')
    end
end
