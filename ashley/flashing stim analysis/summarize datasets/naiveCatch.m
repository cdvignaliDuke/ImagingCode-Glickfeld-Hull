cdirs_all = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)   
        cdirs_all = unique([cdirs_all mouse(imouse).expt(iexp).info.cDirs]);
    end
end

rH = cell(size(cdirs_all));
rCR = cell(size(cdirs_all));
ori_ind = cell(1,4);
rH_all = [];
rCR_all = [];
r_all = [];
rH_cyc = cell(1,2);
rCR_cyc = cell(1,2);
short = 1:4;
% med = 5:6;
long = 5:10;

ncells_ind = 0;

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)        
        if mouse(imouse).expt(iexp).info.isCatch
            cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
%             cell_ind = intersect(mouse(imouse).expt(iexp).cells(13).ind,cell_ind);
            cdirs = mouse(imouse).expt(iexp).info.cDirs;
            cdirs_ind = find(ismember(cdirs_all,cdirs));
            for iori = 1:4
                ori_temp = find(ismember(cell_ind, mouse(imouse).expt(iexp).cells(iori+1).ind))+ncells_ind;
                ori_ind{iori} = cat(1,ori_ind{iori}, ori_temp);
            end
            ncells_ind = ncells_ind+length(cell_ind);
            
            cycH = double(unique(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).tcyc{end}));
            cycCR = double(unique(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).tcyc{end}));
            if min(cycH) > min(cycCR)
                minCyc = min(cycH);
            elseif min(cycH) < min(cycCR)
                minCyc = min(cycCR);
            end
            if max(cycH) > max(cycCR)
                maxCyc = max(cycCR);
            elseif max(cycH) < max(cycCR)
                maxCyc = max(cycH);
            end
            H_temp = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp;
            CR_temp = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp;
            ncycH_temp = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).tcyc;
            ncycCR_temp = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).tcyc;
            rH(cdirs_ind) = cellfun(@(x,y) cat(2,x,y), rH(cdirs_ind), cellfun(@(z,a) mean(z(:,cell_ind,a>minCyc & a<maxCyc),3),H_temp,ncycH_temp,'unif',false),'unif',false);
            rCR(cdirs_ind) = cellfun(@(x,y) cat(2,x,y), rCR(cdirs_ind), cellfun(@(z,a) mean(z(:,cell_ind,a>minCyc & a<maxCyc),3),CR_temp,ncycCR_temp,'unif',false),'unif',false);
            
            H_temp = cat(3,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{:});
            CR_temp = cat(3,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{:});
            ncycH_temp = cat(2,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).tcyc{:});
            ncycCR_temp = cat(2,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).tcyc{:});
            rH_all = cat(2,rH_all,mean(H_temp(:,cell_ind,ncycH_temp>minCyc & ncycH_temp<maxCyc),3));
            rCR_all = cat(2,rCR_all,mean(CR_temp(:,cell_ind,ncycCR_temp>minCyc & ncycCR_temp<maxCyc),3));
            r_all = cat(2,r_all,mean(cat(3,H_temp(:,cell_ind,ncycH_temp>minCyc & ncycH_temp<maxCyc),CR_temp(:,cell_ind,ncycCR_temp>minCyc & ncycCR_temp<maxCyc)),3));
            
            cycs = unique([cycH cycCR]);
            tCycH = double(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).tcyc{end});
            tCycCR = double(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).tcyc{end});
            
            rH_cyc{1} = cat(2,rH_cyc{1},mean(H_temp(:,cell_ind,find(ismember(tCycH,short))),3));
%             rH_cyc{2} = cat(2,rH_cyc{2},mean(H_temp(:,cell_ind,find(ismember(tCycH,med))),3)); 
            rH_cyc{2} = cat(2,rH_cyc{2},mean(H_temp(:,cell_ind,find(ismember(tCycH,long))),3)); 
            
            rCR_cyc{1} = cat(2,rCR_cyc{1},mean(CR_temp(:,cell_ind,find(ismember(tCycCR,short))),3));
%             rCR_cyc{2} = cat(2,rCR_cyc{2},mean(CR_temp(:,cell_ind,find(ismember(tCycCR,med))),3)); 
            rCR_cyc{2} = cat(2,rCR_cyc{2},mean(CR_temp(:,cell_ind,find(ismember(tCycCR,long))),3)); 
        end
    end
end

cyc_name = {'short'; 'long'};

rH_tc = cellfun(@(x) mean(x,2), rH,'unif',false);
rH_tc = cellfun(@(x,y) x-y, rH_tc, cellfun(@(x) mean(x(pre_win,:),1),rH_tc,'unif',false),'unif',false);
rCR_tc = cellfun(@(x) mean(x,2),rCR,'unif',false);
rCR_tc = cellfun(@(x,y) x-y, rCR_tc, cellfun(@(x) mean(x(pre_win,:),1),rCR_tc,'unif',false),'unif',false);

rH_resp = cellfun(@(x,y) y-x,cellfun(@(z) mean(z(pre_win,:),1),rH,'unif',false),cellfun(@(z) mean(z(trans_win,:),1),rH,'unif',false),'unif',false);
rCR_resp = cellfun(@(x,y) y-x,cellfun(@(z) mean(z(pre_win,:),1),rCR,'unif',false),cellfun(@(z) mean(z(trans_win,:),1),rCR,'unif',false),'unif',false);

errH_tc = cellfun(@(x) std(x(:,cell_ind),[],2)/sqrt(size(x,2)),rH,'unif',false);
errCR_tc = cellfun(@(x) std(x(:,cell_ind),[],2)/sqrt(size(x,2)),rCR,'unif',false);

%% plot tc of each dir resp
HvsCR_TC = figure;
suptitle('all val(black) and inval(cyan) targets')
plots = [1:2:length(cdirs_all)*2];
for iplot = 1:length(cdirs_all)
subplot(length(cdirs_all),2,plots(iplot))
P1 = shadedErrorBar(ttMs,rH_tc{iplot},errH_tc{iplot},'k');
hold on
P2 = shadedErrorBar(ttMs,rCR_tc{iplot},errCR_tc{iplot},'c');
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title(num2str(chop(cdirs(iplot),2)))
xlabel('ms')
ylabel('dF/F')
end

%plot cdf of response to target
plots = [2:2:length(cdirs_all)*2];
for iplot = 1:length(cdirs_all)
subplot(length(cdirs_all),2,plots(iplot))
P1 = cdfplot(rH_resp{iplot})
P1.Color = 'k';
hold on
P2 = cdfplot(rCR_resp{iplot})
P2.Color = 'c';
hold on
[h p] = kstest2(rH_resp{iplot},rCR_resp{iplot});
xlim([-0.05 0.05])
title([num2str(chop(cdirs(iplot),2)) '; p= ' num2str(p)])
end

print([fnout 'catch_align_respbytargettype_TCandCDF' datasetStr '.pdf'], '-dpdf')
%% all cells plots
HvsCR_TC = figure;
suptitle('all val(black) and inval(cyan) targets')

subplot(1,2,1)
offset = mean(mean(rH_all(pre_win,:),2),1);
P1 = shadedErrorBar(ttMs,mean(rH_all,2)-offset,std(rH_all,[],2)/sqrt(size(rH_all,2)),'k');
hold on
offset = mean(mean(rCR_all(pre_win,:),2),1);
P2 = shadedErrorBar(ttMs,mean(rCR_all,2)-offset,std(rCR_all,[],2)/sqrt(size(rCR_all,2)),'c');
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title('time-course')
xlabel('ms')
ylabel('dF/F')

subplot(1,2,2)
P1_data = mean(rH_all(trans_win,:),1)-mean(rH_all(pre_win,:),1);
P1 = cdfplot(P1_data);
P1.Color = 'k';
hold on
P2_data = mean(rCR_all(trans_win,:),1)-mean(rCR_all(pre_win,:),1);
P2 = cdfplot(P2_data);
P2.Color = 'c';
hold on
[h p] = kstest2(P1_data,P2_data);
xlim([-0.05 0.05])
title(['kstest, p= ' num2str(p)])

print([fnout 'catch_align_all_TCandCDF' datasetStr '.pdf'], '-dpdf')
% plot by tuning pref
ori_HvsCR = figure;
ori_name = {'0';'45';'90';'135'};
for iori = 1:4
subplot(2,2,iori)
cells = ori_ind{iori};
offset = mean(mean(rH_all(pre_win,cells),2),1);
P1 = shadedErrorBar(ttMs,mean(rH_all(:,cells),2)-offset,std(rH_all(:,cells),[],2)/sqrt(size(rH_all(:,cells),2)),'k');
hold on
offset = mean(mean(rCR_all(pre_win,cells),2),1);
P2 = shadedErrorBar(ttMs,mean(rCR_all(:,cells),2)-offset,std(rCR_all(:,cells),[],2)/sqrt(size(rCR_all(:,cells),2)),'c');
xlim([-10 20]/(cycTime/cycTimeMs))
% ylim([-0.01 0.05])
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title([ori_name{iori} '-slctv;n=' num2str(length(ori_ind{iori}))])
xlabel('ms')
ylabel('dF/F')
end

print([fnout 'catch_align_oriSlctv_TCandCDF' datasetStr '.pdf'], '-dpdf')
%% plot tc of selective cells
cell_ind_pref = cell(1,4);
L = 0;
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
            cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
            cell_ind = intersect(mouse(imouse).expt(iexp).cells(13).ind,cell_ind);   
        for i = 1:4
        cell_ind_pref{i} = [cell_ind_pref{i} find(ismember(cell_ind,mouse(imouse).expt(iexp).cells(i+1).ind))'+L];
        end
            L = length(cell_ind)+L;
    end
end


dirHvsCR = figure;
for i = 1:4
    rH_tc_preftemp = cellfun(@(x) mean(x(:,cell_ind_pref{i}),2),rH,'unif',false);
    rCR_tc_preftemp = cellfun(@(x) mean(x(:,cell_ind_pref{i}),2),rCR,'unif',false);
    rH_prewin = cellfun(@(x) mean(x(pre_win,:)),rH_tc_preftemp,'unif',false);
    rCR_prewin = cellfun(@(x) mean(x(pre_win,:)),rCR_tc_preftemp,'unif',false);
    
%     errH_tc_preftemp = cellfun(@(x) std(x(:,cell_ind_pref{i}),[],2)/sqrt(size(x(:,cell_ind_pref{i}),2)),cellfun(@(x) mean(x,3),rH,'unif',false),'unif',false);
%     errCR_tc_preftemp = cellfun(@(x) std(x(:,cell_ind_pref{i}),[],2)/sqrt(size(x(:,cell_ind_pref{i}),2)),cellfun(@(x) mean(x,3),rCR,'unif',false),'unif',false);

    subplot(2,4,i)
%     P1 = shadedErrorBar(ttMs,rH_tc_preftemp{1},errH_tc_preftemp{1},'k');
%     hold on
%     P2 = shadedErrorBar(ttMs,rCR_tc_preftemp{1},errCR_tc_preftemp{1},'c');
%     hold on
    P1 = plot(ttMs,rH_tc_preftemp{1}-rH_prewin{1},'k');
    hold on
    P2 = plot(ttMs,rCR_tc_preftemp{1}-rCR_prewin{1},'c');
    hold on
    title([mouse(imouse).expt(iexp).cells(i+1).name '-pref; ' num2str(chop(cdirs(1),2)) '-stim'])
    
    subplot(2,4,i+4)
%     P1 = shadedErrorBar(ttMs,rH_tc_preftemp{2},errH_tc_preftemp{2},'k');
%     hold on
%     P2 = shadedErrorBar(ttMs,rCR_tc_preftemp{2},errCR_tc_preftemp{2},'c');
%     hold on
    P1 = plot(ttMs,rH_tc_preftemp{2}-rH_prewin{2},'k');
    hold on
    P2 = plot(ttMs,rCR_tc_preftemp{2}-rCR_prewin{2},'c');
    hold on
    title([mouse(imouse).expt(iexp).cells(i+1).name '-pref; ' num2str(chop(cdirs(2),2)) '-stim'])
        
end 

for i = 1:8
subplot(2,4,i)
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.08])
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
end

%% plot by trial length

toprow = 1:2;
botrow = 3:4;

%short vs long
figure;

val = cellfun(@(x) mean(x,2)-mean(mean(x(pre_win,:),2),1),  rH_cyc,'unif',false);
inv = cellfun(@(x) mean(x,2)-mean(mean(x(pre_win,:),2),1),  rCR_cyc,'unif',false);
val_err = cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),  rH_cyc,'unif',false);
inv_err = cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),  rCR_cyc,'unif',false);

subplot(2,2,1)
P1 = shadedErrorBar(ttMs,val{1},val_err{1},'k',1);
P1.mainLine.LineWidth = 1.5;
P_leg(1) = P1.mainLine;
hold on
P2 = shadedErrorBar(ttMs,val{2},val_err{2},'k',1);
P2.mainLine.Color = [0.5 0.5 0.5];
P2.mainLine.LineWidth = 1.5;
P2.patch.FaceColor = P2.patch.FaceColor*0.5;
P_leg(2) = P2.mainLine;
legend(P_leg,cyc_name)
title('val')
subplot(2,2,2)
P1 = shadedErrorBar(ttMs,inv{1},inv_err{1},'c');
P1.mainLine.LineWidth = 1.5;
P_leg(1) = P1.mainLine;
hold on
P2 = shadedErrorBar(ttMs,inv{2},inv_err{2},'c');
P2.mainLine.Color = 0.5*P2.mainLine.Color;
P2.mainLine.LineWidth = 1.5;
P2.patch.FaceColor = P2.patch.FaceColor*0.5;
P_leg(2) = P2.mainLine;
legend(P_leg,cyc_name)
title('inv')


val = cellfun(@(x) mean(x(trans_win,:),1)-mean(x(pre_win,:),1),  rH_cyc,'unif',false);
inv = cellfun(@(x) mean(x(trans_win,:),1)-mean(x(pre_win,:),1),  rCR_cyc,'unif',false);

subplot(2,2,3)
P1 = cdfplot(val{1})
P1.Color = 'k';
hold on
P2 = cdfplot(val{2})
P2.Color = [0.5 0.5 0.5];
[h p(1)] = kstest2(val{1},val{2});

subplot(2,2,4)
P1 = cdfplot(inv{1})
P1.Color = 'c';
hold on
P2 = cdfplot(inv{2})
P2.Color = P1.Color*.5;
[h p(2)] = kstest2(inv{1},inv{2});

for iplot = 1:2
    subplot(2,2,toprow(iplot))
    xlim([-10 20]/(cycTime/cycTimeMs))
    ylim([-0.01 0.05])
    vline(baseStimFrames/(cycTime/cycTimeMs),':k')
    hold on
    vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
    hold on
    vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
    xlabel('ms')
    ylabel('dF/F')
    
    subplot(2,2,botrow(iplot))
    xlim([-0.05 0.05])
    title(['; p= ' num2str(p(iplot))])
end

print([fnout 'catch_align_longVsShort_lengthcompare_TCandCDF' datasetStr '.pdf'], '-dpdf')

%val vs inv
tarRespbyTrLength = figure;
for iplot = 1:2
    val = rH_cyc{iplot};
    inv = rCR_cyc{iplot};
    
    subplot(2,2,toprow(iplot))
    P1 = shadedErrorBar(ttMs,mean(val,2)-mean(mean(val(pre_win,:),2),1),std(val,[],2)/sqrt(size(val,2)),'k');
    hold on
    P2 = shadedErrorBar(ttMs,mean(inv,2)-mean(mean(inv(pre_win,:),2),1),std(inv,[],2)/sqrt(size(inv,2)),'c');
    xlim([-10 20]/(cycTime/cycTimeMs))
    ylim([-0.01 0.05])
    vline(baseStimFrames/(cycTime/cycTimeMs),':k')
    hold on
    vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
    hold on
    vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
    title(cyc_name{iplot})
    xlabel('ms')
    ylabel('dF/F')
    
    subplot(2,2,botrow(iplot))
    P1 = cdfplot(mean(val(trans_win,:),1)-mean(val(pre_win,:),1));
    P1.Color = 'k';
    hold on
    P2 = cdfplot(mean(inv(trans_win,:),1)-mean(inv(pre_win,:),1))
    P2.Color = 'c';
    hold on
    [h p] = kstest2(mean(val(trans_win,:),1)-mean(val(pre_win,:),1),mean(inv(trans_win,:),1)-mean(inv(pre_win,:),1));
    xlim([-0.05 0.05])
    title(['; p= ' num2str(p)])
end

print([fnout 'catch_align_longVsShort_typecompare_TCandCDF' datasetStr '.pdf'], '-dpdf')
