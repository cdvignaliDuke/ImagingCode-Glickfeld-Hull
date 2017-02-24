% make bins across trial lengths
trL_val = [];

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)        
        if mouse(imouse).expt(iexp).info.isCatch  
            cycTimeMs_temp = mouse(imouse).expt(iexp).info.cyc_time_ms;
            cycTime_temp = mouse(imouse).expt(iexp).info.cyc_time;
            
            cT = cat(2,cat(2,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).tcyc{hInd}),cat(2,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).tcyc{mInd}));
            msT = cT*cycTimeMs_temp;
            
            trL_val = cat(2,trL_val,msT);            
        end
    end
end

trL = unique(trL_val);
            
bin_edges = double([0 2000 3000 max(trL)]);
[n bins] = histc(double(trL_val),bin_edges);
nbins = length(unique(bins));
% create variables with all valid and all invalid trials trials (from catch experiments only) and matched size variables with
% trial length
cellsCatchAlign = 13;

tc_val = cell(1,nbins);
trL_val_binned = cell(1,nbins);
tc_inv = cell(1,nbins);
trL_inv_binned = cell(1,nbins);

resp_eachTr_val = [];
trL_eachTr_val = [];
resp_eachTr_inv = [];
trL_eachTr_inv = [];

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)        
        if mouse(imouse).expt(iexp).info.isCatch            
            cell_ind = mouse(imouse).expt(iexp).cells(cellsCatchAlign).ind;
            cycTimeMs_temp = mouse(imouse).expt(iexp).info.cyc_time_ms;
            cycTime_temp = mouse(imouse).expt(iexp).info.cyc_time;
            
            %valid trials
            hInd = find(cell2mat(cellfun(@(x) ~isempty(x),mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp,'unif',false)));
            mInd = find(cell2mat(cellfun(@(x) ~isempty(x),mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp,'unif',false)));
            
            rT = cat(3,cat(3,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{hInd}),cat(3,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp{mInd}));
            msT = cycTimeMs_temp*(cat(2,cat(2,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).tcyc{hInd}),cat(2,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).tcyc{mInd})));
            
            resp_eachTr_val = cat(1,resp_eachTr_val,squeeze(bsxfun(@minus,mean(mean(rT(trans_win,cell_ind,:),2),1),mean(mean(rT(pre_win,cell_ind,:),2),1))));
            trL_eachTr_val = cat(2,trL_eachTr_val,msT);
            
            [n bins] = histc(double(msT),bin_edges);
            bin_ind = find(n ~= 0);
            for ibin = 1:length(unique(bins))
                trInd = find(bins == ibin);
                tc_val{bin_ind(ibin)} = cat(2,tc_val{bin_ind(ibin)},mean(rT(:,cell_ind,trInd),3));
                trL_val_binned{bin_ind(ibin)} = cat(2,trL_val_binned{bin_ind(ibin)},msT(trInd));                
            end
            
            %invalid trials
            faInd = find(cell2mat(cellfun(@(x) ~isempty(x),mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp,'unif',false)));
            crInd = find(cell2mat(cellfun(@(x) ~isempty(x),mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp,'unif',false)));
            
            rT = cat(3,cat(3,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{faInd}),cat(3,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{crInd}));
            msT = cycTimeMs_temp*(cat(2,cat(2,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).tcyc{faInd}),cat(2,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).tcyc{crInd})));

            resp_eachTr_inv = cat(1,resp_eachTr_inv,squeeze(bsxfun(@minus,mean(mean(rT(trans_win,cell_ind,:),2),1),mean(mean(rT(pre_win,cell_ind,:),2),1))));
            trL_eachTr_inv = cat(2,trL_eachTr_inv,msT);
            
            [n bins] = histc(double(msT),bin_edges);
            bin_ind = find(n ~= 0);
            for ibin = 1:length(unique(bins))
                trInd = find(bins == ibin);
                tc_inv{bin_ind(ibin)} = cat(2,tc_inv{bin_ind(ibin)},mean(rT(:,cell_ind,trInd),3));
                trL_inv_binned{bin_ind(ibin)} = cat(2,trL_inv_binned{bin_ind(ibin)},msT(trInd));                
            end
            
        end
    end
end
trL_val_mean = chop(cell2mat(cellfun(@mean,trL_val_binned,'unif',0)),3);
tc_val_mean = cellfun(@(x) bsxfun(@minus,mean(x,2),mean(mean(x(pre_win,:),2),1)),tc_val,'unif',0);
tc_val_ste = cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),tc_val,'unif',0);
resp_val_mean = cellfun(@(x) bsxfun(@minus,mean(x(trans_win,:),1),mean(x(pre_win,:),1)),tc_val,'unif',0);

tc_inv_emptybins = find(cell2mat(cellfun(@(x) lt(length(x),5),trL_inv_binned,'unif',0)));
nonempties = setdiff(1:nbins,tc_inv_emptybins);
trL_inv_mean = chop(cell2mat(cellfun(@mean,trL_inv_binned(nonempties),'unif',0)),3);
tc_inv_mean = cellfun(@(x) bsxfun(@minus,mean(x,2),mean(mean(x(pre_win,:),2),1)),tc_inv(nonempties),'unif',0);
tc_inv_ste = cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),tc_inv(nonempties),'unif',0);
resp_inv_mean = cellfun(@(x) bsxfun(@minus,mean(x(trans_win,:),1),mean(x(pre_win,:),1)),tc_inv(nonempties),'unif',0);


% plot tc,cdfs, and scatters of responses for each time-bin
colorsVal = brewermap(nbins+1,'Greys');
colorsVal = colorsVal(2:end,:);

colorsInv = brewermap(length(nonempties)+1,'Blues');
colorsInv = colorsInv(2:end,:);

%time-course
valTimeBinnedFig = figure;
subplot(2,3,1)
for ibin = 1:nbins
h = shadedErrorBar(ttMs,tc_val_mean{ibin},tc_val_ste{ibin},'k');
hold on
h.mainLine.Color = colorsVal(ibin,:);
h.mainLine.LineWidth = 3;
leg(ibin) = h.mainLine;
end
xlim([-15 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
axis square
legend(leg,strread(num2str(bin_edges),'%s'),'Location','northwest')
title('valid')
xlabel('time (ms)')
ylabel('dF/F')

subplot(2,3,4)
for ibin = 1:length(nonempties)
h = shadedErrorBar(ttMs,tc_inv_mean{nonempties(ibin)},tc_inv_ste{nonempties(ibin)},'k');
hold on
h.mainLine.Color = colorsInv(ibin,:);
h.mainLine.LineWidth = 3;
leg(ibin) = h.mainLine;
end
xlim([-15 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
axis square
legend(leg,strread(num2str(bin_edges(nonempties)),'%s'),'Location','northwest')
title('invalid')
xlabel('time (ms)')
ylabel('dF/F')

%cdf
subplot(2,3,2)
for ibin = 1:nbins
h = cdfplot(resp_val_mean{ibin});
h.Color = colorsVal(ibin,:);
hold on
leg2(ibin) = h;
end
xlim([-0.05 0.05])
axis square
legend(leg2,strread(num2str(trL_val_mean),'%s'),'location','northwest')
xlabel('dF/F')
ylabel('fraction of cells')

subplot(2,3,5)
for ibin = 1:length(nonempties)
h = cdfplot(resp_inv_mean{nonempties(ibin)});
h.Color = colorsInv(ibin,:);
hold on
leg2(ibin) = h;
end
xlim([-0.05 0.05])
axis square
legend(leg2,strread(num2str(trL_inv_mean(nonempties)),'%s'),'location','northwest')
xlabel('dF/F')
ylabel('fraction of cells')

%scatter
subplot(2,3,6)
h = scatter(trL_eachTr_val,resp_eachTr_val, 'o','k');
h.MarkerFaceColor = [0 0 0];
hold on
h = scatter(trL_eachTr_inv,resp_eachTr_inv, 'o','c');
h.MarkerFaceColor = [0 1 1];
hold on
axis square
xlabel('trial length (ms)')
ylabel('resp (dF/F)')

figure(valTimeBinnedFig);
print([fnout '_trialLength_val-inv-all'],'-dpdf','-fillpage')