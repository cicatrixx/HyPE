function CompileMinMaxByRcpGroups(myhistdatalumped, myfutdatalumped, oxlsname)
%% Compile min, max and mean for hist and fut TS grouped by RCP and corners
% datalumped are 2 param structs w first col for technical pf and second col for
% sustainable pf
% each param contains a 3D matrix for pf x dataseries x qf

run('myVarNames_Fut.m')
histidx=1;
selnewpots=[1 4];

for pot=1:2
    colstats_nf=[];
    colstats_ff=[];
    colstats_hist_nf=[];
    colstats_hist_ff=[];

    myhistdata=myhistdatalumped{pot};
    myfutdata=myfutdatalumped{pot};
    
    %hist pf under histQ - 1 row added
    curdata=myhistdata(histidx,:,1);
    colstats_hist_nf=[colstats_hist_nf;  mean(curdata) min(curdata) max(curdata)];

    for ssel=1:3
        %select two futures TSs
        qselnf=sspgrps{ssel};
        qselff=sspgrps{ssel+3};

        % add  hist pf under rcpQs NF and FF tp top of pile
        % there are 4 corner Qs for NF - take mean of them  - 1 row added
        %nf
        curdata=squeeze(myfutdata(histidx,:,qselnf));
        colstats_hist_nf=[ colstats_hist_nf; mean(curdata,'all')' min(curdata,[],'all')' max(curdata,[],'all')'];
        %ff
        curdata=squeeze(myfutdata(histidx,:,qselff));
        colstats_hist_ff=[ colstats_hist_ff; mean(curdata,'all')' min(curdata,[],'all')' max(curdata,[],'all')' ];

        %  4 corner pf under rcpQs for NF and FF - 4 rows added - one for
        %  each corner
        for cr=1:4
            %nf
            curdata=squeeze(myfutdata(qselnf(cr),:,qselnf));
            colstats_nf=[colstats_nf;  mean(curdata,'all')' min(curdata,[],'all')' max(curdata,[],'all')'];

            %ff
            curdata=squeeze(myfutdata(qselff(cr),:,qselff));
            colstats_ff=[colstats_ff;  mean(curdata,'all')' min(curdata,[],'all')' max(curdata,[],'all')'];
        end
    end
    gr_temporalstats_futbyrcp{pot}=[colstats_nf;colstats_ff];
    gr_temporalstats_histbyrcp{pot}=[colstats_hist_nf;colstats_hist_ff];
end
disp("Eval min max growth rate")

%% Convert to table, eval mean range for fut TS grouped into 6 (3RCP x 2TF)
for pot=1:2
    curpotname=newpots4{selnewpots(pot)};
    tblfut{pot}=array2table(gr_temporalstats_futbyrcp{pot},"VariableNames",strcat(curpotname,"_",["Mean" "Min" "Max"]));
    tblfut{pot}.pfnames=tfhistrcpcornernames(2:end)';
        meanrange{pot}=groupsummary(tblfut{pot}{:,1},repelem(1:6,1,4)',{'min' 'max'});
    tblhist{pot}=array2table(gr_temporalstats_histbyrcp{pot},"VariableNames",strcat(curpotname,"_",["Mean" "Min" "Max"]));
    tblhist{pot}.Qnames=tfhistrcpnames';
end
tmptbl=table(meanrange{1},meanrange{2},'VariableNames',newpots4(selnewpots(1:2))); %"mean-min" "mean-max"
tmptbl.scnames=tfrcpnames';

%% Write growthrates to excel
oxlsfile=fullfile(rootoffut,oxlsname);

% Fut PF under Fut Qs
writetable([tblfut{1}(:,[4 1:3]) tblfut{2}(:,1:3)],oxlsfile,'Sheet',"fut")
% Hist PF under Hist+Fut Qs
writetable([tblhist{1}(:,[4 1:3]) tblhist{2}(:,1:3)],oxlsfile,'Sheet',"hist")
% Range in Fut PF_FutQ
writetable(tmptbl,oxlsfile,'Sheet',"fut_mean_range")
disp("Written to excel")

end