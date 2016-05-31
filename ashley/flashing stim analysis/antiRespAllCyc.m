
cycs = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        cycstemp = 1:size(mouse(imouse).expt(iexp).align(1).av(1).outcome(1).cmlvCycResp,2);
        cycs = unique([cycs cycstemp]);
    end
end
ncyc = length(cycs);

c = ceil(sqrt(ncyc+1));
if (c^2)-c > ncyc+1
   c2 = c-1;
else
    c2 = c;
end
%% cumulative trials upto number of trials
% i = 1;
resp_vis_cmlvCyc = cell(1,length(cycs));
resp_aud_cmlvCyc = cell(1,length(cycs));
resp_all_cmlvCyc = cell(1,length(cycs));
% std_vis = [];
% std_aud = [];


for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(1).ind,cell_ind);
        
        for icyc = 1:length(cycs)
        if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp,2) >= icyc
            resp_vis_cmlvCyc{icyc} = cat(2,resp_vis_cmlvCyc{icyc},mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp{icyc}(:,cell_ind,:),3));
            resp_aud_cmlvCyc{icyc} = cat(2,resp_aud_cmlvCyc{icyc},mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvCycResp{icyc}(:,cell_ind,:),3));
            resp_all_cmlvCyc{icyc} = cat(2,resp_all_cmlvCyc{icyc},mean(cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp{icyc}(:,cell_ind,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvCycResp{icyc}(:,cell_ind,:)),3));
        end
        end
            
%        i = i+1;
        end
end

base_win = [pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames];
resp_win = [trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames];

resp_vis_cmlvCyc_resp = cell(1,length(cycs));
resp_aud_cmlvCyc_resp= cell(1,length(cycs));
for icyc = 1:length(cycs)
    r = (trans_win(1):trans_win(end))+((icyc-1)*cycTime);
    b = (pre_win(1):pre_win(end))+((icyc-1)*cycTime);
   resp_vis_cmlvCyc_resp{icyc} = bsxfun(@minus, mean(resp_vis_cmlvCyc{icyc}(r,:),1), mean(resp_vis_cmlvCyc{icyc}(b,:),1)); 
   resp_aud_cmlvCyc_resp{icyc} = bsxfun(@minus, mean(resp_aud_cmlvCyc{icyc}(r,:),1), mean(resp_aud_cmlvCyc{icyc}(b,:),1));     
end
nCells = size(resp_vis_cmlvCyc_resp{1},2);

pcmlvCyc = zeros(1,length(cycs));
for icyc = 1:length(cycs)
[h, pcmlvCyc(icyc)] = ttest(resp_vis_cmlvCyc_resp{icyc}, resp_aud_cmlvCyc_resp{icyc},'alpha', 0.05/nCells);%, 'dim', 2, 'tail', 'left', 'alpha', 0.05/nCells);
end

resp_vis_cmlvCyc_mean = cellfun(@(x) mean(x,2),resp_vis_cmlvCyc,'unif',false);
resp_aud_cmlvCyc_mean = cellfun(@(x) mean(x,2),resp_aud_cmlvCyc,'unif',false);

ttCycMs = (-pre_event_frames:cycs(end)*cycTime)/(cycTime/cycTimeMs);
respAvsVFig_cmlvCycs = figure;
suptitle({[titleStr '; alpha = ' num2str(0.05/nCells)]})
for icyc = 1:length(cycs)
figure(respAvsVFig_cmlvCycs);
subplot(c,c2,icyc)
plot(ttCycMs(1:size(resp_vis_cmlvCyc{icyc},1)), resp_vis_cmlvCyc_mean{icyc}, 'g')
hold on
plot(ttCycMs(1:size(resp_vis_cmlvCyc{icyc},1)), resp_aud_cmlvCyc_mean{icyc}, 'k')
vline((resp_win+((icyc-1)*cycTime))/(cycTime/cycTimeMs), '--r')
vline((base_win+((icyc-1)*cycTime))/(cycTime/cycTimeMs), '--k')
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
xlim([-10 cycTime*icyc]/(cycTime/cycTimeMs))
% ylim([-0.01 0.03]);
title([num2str(icyc) ' cycs;All cells; p = ' num2str(pcmlvCyc(icyc))])
end

resp_vis_cmlvCyc_ste = cellfun(@(x) x/sqrt(nCells), cellfun(@(x) std(x,[],2), resp_vis_cmlvCyc,'unif',false), 'unif',false) ;
resp_aud_cmlvCyc_ste = cellfun(@(x) x/sqrt(nCells), cellfun(@(x) std(x,[],2), resp_aud_cmlvCyc,'unif',false), 'unif',false) ;

respAvsVFig_cmlvCycsErr = figure;
suptitle({[titleStr '; alpha = ' num2str(0.05/nCells)]})
for icyc = 1:length(cycs)
figure(respAvsVFig_cmlvCycsErr);
subplot(c,c2,icyc)
shadedErrorBar(ttCycMs(1:size(resp_vis_cmlvCyc{icyc},1)), resp_vis_cmlvCyc_mean{icyc},resp_vis_cmlvCyc_ste{icyc},'g');
hold on
shadedErrorBar(ttCycMs(1:size(resp_aud_cmlvCyc{icyc},1)), resp_aud_cmlvCyc_mean{icyc},resp_aud_cmlvCyc_ste{icyc},'k');

vline((resp_win+((icyc-1)*cycTime)/(cycTime/cycTimeMs)), '--r')
vline((base_win+((icyc-1)*cycTime)/(cycTime/cycTimeMs)), '--k')
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
xlim([-10 cycTime*icyc]/(cycTime/cycTimeMs))
% ylim([-0.01 0.03]);
title([num2str(icyc) ' cycs;All cells; p = ' num2str(pcmlvCyc(icyc))])
end

figure(respAvsVFig_cmlvCycs);
print([fnout 'press_align_TCbycmlvCyc' datasetStr '.pdf'], '-dpdf');
figure(respAvsVFig_cmlvCycsErr);
print([fnout 'press_align_TCbycmlvCycErr' datasetStr '.pdf'], '-dpdf')

%% trials of a specific length
resp_vis_cyc = cell(1,length(cycs));
resp_aud_cyc = cell(1,length(cycs));
resp_all_cyc = cell(1,length(cycs));
nTrials_vis = zeros(1,length(cycs));
nTrials_aud = zeros(1,length(cycs));
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(1).ind,cell_ind);
        
        for icyc = 1:length(cycs)
            if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp,2) >= icyc
        if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp{icyc},3) > 0 
            resp_vis_cyc{icyc} = cat(2,resp_vis_cyc{icyc},mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{icyc}(:,cell_ind,:),3));
            resp_aud_cyc{icyc} = cat(2,resp_aud_cyc{icyc},mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cycResp{icyc}(:,cell_ind,:),3));
            resp_all_cyc{icyc} = cat(2,resp_all_cyc{icyc},mean(cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{icyc}(:,cell_ind,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cycResp{icyc}(:,cell_ind,:)),3));
            nTrials_vis(icyc) = nTrials_vis(icyc)+size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{icyc}(:,cell_ind,:),3);
            nTrials_aud(icyc) = nTrials_aud(icyc)+size(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cycResp{icyc}(:,cell_ind,:),3);
        end
            end
        end
            
%        i = i+1;
        end
end

base_win = [pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames];
resp_win = [trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames];

resp_vis_cyc_resp = cell(1,length(cycs));
resp_aud_cyc_resp= cell(1,length(cycs));
for icyc = 1:length(cycs)
    r = (trans_win(1):trans_win(end))+((icyc-1)*cycTime);
    b = (pre_win(1):pre_win(end))+((icyc-1)*cycTime);
   resp_vis_cyc_resp{icyc} = bsxfun(@minus, mean(resp_vis_cyc{icyc}(r,:),1), mean(resp_vis_cyc{icyc}(b,:),1)); 
   resp_aud_cyc_resp{icyc} = bsxfun(@minus, mean(resp_aud_cyc{icyc}(r,:),1), mean(resp_aud_cyc{icyc}(b,:),1));     
end
nCells = size(resp_vis_cyc_resp{1},2);

pCyc = zeros(1,length(cycs));
for icyc = 1:length(cycs)
[h, pCyc(icyc)] = ttest(resp_vis_cyc_resp{icyc}, resp_aud_cyc_resp{icyc},'alpha', 0.05/nCells);%, 'dim', 2, 'tail', 'left', 'alpha', 0.05/nCells);
end

resp_vis_cyc_mean = cellfun(@(x) nanmean(x,2),resp_vis_cyc,'unif',false);
resp_aud_cyc_mean = cellfun(@(x) nanmean(x,2),resp_aud_cyc,'unif',false);

ttCycMs = (-pre_event_frames:cycs(end)*cycTime)/(cycTime/cycTimeMs);
respAvsVFig_cycs = figure;
suptitle({[titleStr '; alpha = ' num2str(0.05/nCells)]})
for icyc = 1:length(cycs)
figure(respAvsVFig_cycs);
subplot(c,c2,icyc)
plot(ttCycMs(1:size(resp_vis_cyc{icyc},1)), resp_vis_cyc_mean{icyc}, 'g')
hold on
plot(ttCycMs(1:size(resp_vis_cyc{icyc},1)), resp_aud_cyc_mean{icyc}, 'k')
vline((resp_win+((icyc-1)*cycTime))/(cycTime/cycTimeMs), '--r')
vline((base_win+((icyc-1)*cycTime))/(cycTime/cycTimeMs), '--k')
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
xlim([-10 cycTime*icyc]/(cycTime/cycTimeMs))
% ylim([-0.01 0.03]);
title([num2str(icyc) ' cycs;All cells; p = ' num2str(pCyc(icyc))])
end

resp_vis_cyc_ste = cellfun(@(x) x/sqrt(nCells), cellfun(@(x) nanstd(x,[],2), resp_vis_cyc,'unif',false), 'unif',false) ;
resp_aud_cyc_ste = cellfun(@(x) x/sqrt(nCells), cellfun(@(x) nanstd(x,[],2), resp_aud_cyc,'unif',false), 'unif',false) ;

respAvsVFig_cycsErr = figure;
suptitle({[titleStr '; alpha = ' num2str(0.05/nCells)]})
for icyc = 1:length(cycs)
figure(respAvsVFig_cycsErr);
subplot(c,c2,icyc)
shadedErrorBar(ttCycMs(1:size(resp_vis_cyc{icyc},1)), resp_vis_cyc_mean{icyc},resp_vis_cyc_ste{icyc},'g');
hold on
shadedErrorBar(ttCycMs(1:size(resp_aud_cyc{icyc},1)), resp_aud_cyc_mean{icyc},resp_aud_cyc_ste{icyc},'k');

vline((resp_win+((icyc-1)*cycTime))/(cycTime/cycTimeMs), '--r')
vline((base_win+((icyc-1)*cycTime))/(cycTime/cycTimeMs), '--k')
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
xlim([-10 cycTime*icyc]/(cycTime/cycTimeMs))
% ylim([-0.01 0.03]);
title([num2str(icyc) ' cycs;' num2str(nTrials_vis(icyc)) '/' num2str(nTrials_aud(icyc)) ' vis/aud trials; p = ' num2str(pCyc(icyc))])
end

figure(respAvsVFig_cycs);
print([fnout 'press_align_TCbyCyc' datasetStr '.pdf'], '-dpdf');
figure(respAvsVFig_cycsErr);
print([fnout 'press_align_TCbyCycErr' datasetStr '.pdf'], '-dpdf')