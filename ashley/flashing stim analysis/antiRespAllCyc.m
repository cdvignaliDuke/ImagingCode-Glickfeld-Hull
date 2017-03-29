
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
disp(datasetStr)
%%
if strcmp(datasetStr,'_V1')
    endTrRespWinMs = 300;
    nCycShort = 1:5;
    nCycLong = nCycShort(end)+1:length(cycs);
nTrials_vis = zeros(1,length(cycs));
nTrials_aud = zeros(1,length(cycs));
endTrResp_vis = cell(1,length(cycs));
endTrResp_aud = cell(1,length(cycs));
endTrResp_vis_short = [];
endTrResp_aud_short = [];
endTrResp_vis_long = [];
endTrResp_aud_long = [];
tc_vis_long = [];
tc_aud_long = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(respCellsInd).ind,cell_ind);
        endTrialRespIntervalFrames = endTrRespWinMs*(mouse(imouse).expt(iexp).info.cyc_time/mouse(imouse).expt(iexp).info.cyc_time_ms);
            
        for icyc = 1:length(cycs)
            if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp,2) >= icyc
        if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{icyc},3) > 0 
            nTrials_vis(icyc) = nTrials_vis(icyc)+size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{icyc}(:,cell_ind,:),3);
            nTrials_aud(icyc) = nTrials_aud(icyc)+size(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cycResp{icyc}(:,cell_ind,:),3);
            endTrResp_vis{icyc} = cat(2,endTrResp_vis{icyc},mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{icyc}(end-endTrialRespIntervalFrames:end,cell_ind,:),3),1));
            endTrResp_aud{icyc} = cat(2,endTrResp_aud{icyc},mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cycResp{icyc}(end-endTrialRespIntervalFrames:end,cell_ind,:),3),1));            
%             if sum(ismember(nCycShort,icyc))
%                 endTrResp_vis_short(icyc) = 
%             end
        end
            end
        end
        
        vis_temp = [];
        aud_temp = [];
        for icyc = 1:length(nCycShort)
            if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp,2) >= nCycShort(icyc)
        if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{nCycShort(icyc)},3) > 0 
            vis_temp = cat(3,vis_temp,mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{nCycShort(icyc)}(end-endTrialRespIntervalFrames:end,cell_ind,:),1));            
            aud_temp = cat(3,vis_temp,mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{nCycShort(icyc)}(end-endTrialRespIntervalFrames:end,cell_ind,:),1));
        end
            end
        end
        endTrResp_vis_short = cat(2,endTrResp_vis_short, mean(vis_temp,3));
        endTrResp_aud_short = cat(2,endTrResp_aud_short,mean(aud_temp,3));
        vis_temp = [];
        aud_temp = [];
        tc_v_temp = [];
        tc_a_temp = [];
        for icyc = 1:length(nCycLong)
            if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp,2) >= nCycLong(icyc)
        if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{nCycLong(icyc)},3) > 0 
            if isempty(vis_temp)
                tc_v_temp = mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,cell_ind,:),3);
                tc_a_temp = mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,cell_ind,:),3);
            end
            vis_temp = cat(3,vis_temp,mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{nCycLong(icyc)}(end-endTrialRespIntervalFrames:end,cell_ind,:),1));            
            aud_temp = cat(3,vis_temp,mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{nCycLong(icyc)}(end-endTrialRespIntervalFrames:end,cell_ind,:),1));
        end
            end
        end
        endTrResp_vis_long = cat(2,endTrResp_vis_long, mean(vis_temp,3));
        endTrResp_aud_long = cat(2,endTrResp_aud_long,mean(aud_temp,3));
        tc_vis_long = cat(2,tc_vis_long,tc_v_temp);
        tc_aud_long = cat(2,tc_aud_long,tc_a_temp);
        
    end
end
respByCycScatter = figure;
suptitle([titleStr '-resp from 100ms after last base stim to end'])
for icyc = 1:length(cycs)
    if ~isempty(endTrResp_vis{icyc})
        subplot(c,c2,icyc)
        scatter(endTrResp_vis{icyc},endTrResp_aud{icyc},50,'k.');
        hold on
        plot(-20:1:20,-20:1:20,'k--')
        ax = [min([endTrResp_vis{icyc} endTrResp_aud{icyc}]) max([endTrResp_vis{icyc} endTrResp_aud{icyc}])];
        xlim(ax)
        ylim(ax)
        xlabel('visual')
        ylabel('auditory')
        axis square
        title([num2str(icyc) ' cycs;' num2str(nTrials_vis(icyc)) '/' num2str(nTrials_aud(icyc))])
    end
end
figure(respByCycScatter);
% print([fnout 'press_align_EndTrialsRespbyCyc' datasetStr '.pdf'], '-dpdf');

endTrRespLongVsShort = figure;
suptitle([titleStr '- resp 300ms before target'])
subplot(1,2,1)
scatter(endTrResp_vis_short,endTrResp_aud_short,50,'k.');
hold on
plot(-20:1:20,-20:1:20,'k--')
% ax = [min([endTrResp_vis_long endTrResp_aud_long]) max([endTrResp_aud_long endTrResp_vis_long])];
ax = [-0.1 0.3];
xlim(ax)
ylim(ax)
xlabel('visual')
ylabel('auditory')
axis square
title('short trials - <=5cycs')
subplot(1,2,2)
scatter(endTrResp_vis_long,endTrResp_aud_long,50,'k.');
hold on
plot(-20:1:20,-20:1:20,'k--')
xlim(ax)
ylim(ax)
xlabel('visual')
ylabel('auditory')
axis square
title('long trials - >6cycs')
%% analysis by trial length (ms)
maxTrialLength = 0;
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        if length(mouse(imouse).expt(iexp).align(1).av(1).outcome(1).cmlvCycResp)*mouse(imouse).expt(iexp).info.cyc_time_ms > maxTrialLength
        maxTrialLength = length(mouse(imouse).expt(iexp).align(1).av(1).outcome(1).cmlvCycResp)*mouse(imouse).expt(iexp).info.cyc_time_ms;
        end
    end
end
time_bins = 500:500:maxTrialLength;
for i = 1:5
endTrTimeBin_v(i).oriResp = cell(1,length(time_bins));
endTrTimeBin_a(i).oriResp = cell(1,length(time_bins)); 
end
        
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(respCellsInd).ind,cell_ind);
        time_ms = (-mouse(imouse).expt(iexp).pre_event_frames:max(cell2mat(cellfun(@(x) size(x,1),mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp,'unif',false))))*(mouse(imouse).expt(iexp).info.cyc_time_ms/mouse(imouse).expt(iexp).info.cyc_time);
        sz = cell2mat(cellfun(@(x) size(x,1),mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp,'unif',false))*(mouse(imouse).expt(iexp).info.cyc_time_ms/mouse(imouse).expt(iexp).info.cyc_time);
        trL = sz-(mouse(imouse).expt(iexp).pre_event_frames*(mouse(imouse).expt(iexp).info.cyc_time_ms/mouse(imouse).expt(iexp).info.cyc_time));
            
        for it = 1:length(time_bins)
        if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp,2)*mouse(imouse).expt(iexp).info.cyc_time_ms >= time_bins(it)
            icyc = find(trL > time_bins(it),1);
            if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp{icyc},3) > 0 
                [time_val2 time_ind2] = closest(time_ms,time_bins(it));
                [time_val1 time_ind1] = closest(time_ms,time_bins(it)-endTrRespWinMs);
                endTrTimeBin_v(1).oriResp{it} = cat(2,endTrTimeBin_v(1).oriResp{it},mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp{icyc}(time_ind1:time_ind2,cell_ind,:),3),1));
                endTrTimeBin_a(1).oriResp{it} = cat(2,endTrTimeBin_a(1).oriResp{it},mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvCycResp{icyc}(time_ind1:time_ind2,cell_ind,:),3),1));
                for iori = 1:4
                    cell_ind_ori = intersect(cell_ind,mouse(imouse).expt(iexp).cells(iori+1).ind);
                    if ~isempty(cell_ind_ori)
%                     endTrTimeBin_v(iori+1).oriResp{it} = [];
%                     endTrTimeBin_a(iori+1).oriResp{it} = [];
%                     else
                    endTrTimeBin_v(iori+1).oriResp{it} = cat(2,endTrTimeBin_v(iori+1).oriResp{it},mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp{icyc}(time_ind1:time_ind2,cell_ind_ori,:),3),1));
                    endTrTimeBin_a(iori+1).oriResp{it} = cat(2,endTrTimeBin_a(iori+1).oriResp{it},mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvCycResp{icyc}(time_ind1:time_ind2,cell_ind_ori,:),3),1));
                    end
                end
            end
        end
        end
    end
end

%plot scatters of time bins
tb = ceil(sqrt(ncyc+1));
if (tb^2)-tb > length(time_bins)+1
   tb2 = tb-1;
else
    tb2 = tb;
end
oriCols = [1,0,0; 0,1,0; 0,0,1; 1,0,1];%rgb red, green, blue, magenta
timeBinScatter = figure;
suptitle('end tr mean resp for trial lengths')
for it = 1:length(time_bins)
    subplot(tb,tb2,it)
    allCellsResp = scatter(endTrTimeBin_v(1).oriResp{it},endTrTimeBin_a(1).oriResp{it},100,'k.');
    allCellsResp.CData = [0.7 0.7 0.7];
    hold on
    for iori = 1:4
        oriCellsResp(iori) = scatter(endTrTimeBin_v(iori+1).oriResp{it},endTrTimeBin_a(iori+1).oriResp{it},100,'k.');
        oriCellsResp(iori).CData = oriCols(iori,:);
        hold on
    end
    hold on
    plot(-20:1:20,-20:1:20,'k--')
    axis square
    title([num2str(time_bins(it)) 'ms'])
    xlim([-0.1 0.5])
    ylim([-0.1 0.5])
    xlabel('visual')
    ylabel('auditory')
legend([allCellsResp oriCellsResp],'all','0deg','45deg','90deg','135deg','Location','northeastoutside')
end

%calc modulation index (V-A)/(V+A),plot histogram
for iori = 1:5
    timeBin_mi{iori} = cellfun(@rdivide,cellfun(@minus,endTrTimeBin_v(iori).oriResp,endTrTimeBin_a(iori).oriResp,'unif',false),cellfun(@plus,endTrTimeBin_v(iori).oriResp,endTrTimeBin_a(iori).oriResp,'unif',false),'unif',false);
end
timeBinMIHits = figure;
bin_edges = -20:40/99:20;
suptitle('end tr mod idx for trial lengths')
for iori = 1:5
for it = 1:length(time_bins)
    subplot(tb,tb2,it)
    h = histc(timeBin_mi{iori}{it},bin_edges);
        b = bar(h);
    if iori == 1
        b.FaceColor = [0.7 0.7 0.7];
    else
        b.FaceColor = oriCols(iori-1,:);
    end
    b.EdgeColor = [1 1 1];
    hold on
    
    title([num2str(time_bins(it)) 'ms'])
    xlabel('dF/F')
    ylabel('n cells')
    ylim([0 100])
end
end
%%
else
%% cumulative trials upto number of trials
% i = 1;
resp_vis_cmlvCyc = cell(1,length(cycs));
resp_aud_cmlvCyc = cell(1,length(cycs));
resp_all_cmlvCyc = cell(1,length(cycs));
fano_vis_cmlvCyc = cell(1,length(cycs));
fano_aud_cmlvCyc = cell(1,length(cycs));
fano_all_cmlvCyc = cell(1,length(cycs));
% std_vis = [];
% std_aud = [];

disp('got to resp_vis_cmlvCyc')

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(respCellsInd).ind,cell_ind);
        if mouse(imouse).expt(iexp).info.cyc_time == 11;
            for icyc = 1:length(cycs)
            if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp,2) >= icyc
                resp_vis_cmlvCyc{icyc} = cat(2,resp_vis_cmlvCyc{icyc},mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp{icyc}(:,cell_ind,:),3));
                resp_aud_cmlvCyc{icyc} = cat(2,resp_aud_cmlvCyc{icyc},mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvCycResp{icyc}(:,cell_ind,:),3));
                resp_all_cmlvCyc{icyc} = cat(2,resp_all_cmlvCyc{icyc},mean(cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp{icyc}(:,cell_ind,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvCycResp{icyc}(:,cell_ind,:)),3));
                fano_vis_cmlvCyc{icyc} = cat(2,fano_vis_cmlvCyc{icyc},fano(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp{icyc}(:,cell_ind,:),3));
                fano_aud_cmlvCyc{icyc} = cat(2,fano_aud_cmlvCyc{icyc},fano(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvCycResp{icyc}(:,cell_ind,:),3));
                fano_all_cmlvCyc{icyc} = cat(2,fano_all_cmlvCyc{icyc},fano(cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvCycResp{icyc}(:,cell_ind,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvCycResp{icyc}(:,cell_ind,:)),3));
            end
            end
        end
%        i = i+1;
        end
end
nCells = size(resp_vis_cmlvCyc{1},2);

fano_vis_cmlvCyc_mean = cellfun(@(x) mean(x,2),fano_vis_cmlvCyc,'unif',false);
fano_aud_cmlvCyc_mean = cellfun(@(x) mean(x,2),fano_aud_cmlvCyc,'unif',false);
fano_all_cmlvCyc_mean = cellfun(@(x) mean(x,2),fano_all_cmlvCyc,'unif',false);
fano_vis_cmlvCyc_ste = cellfun(@(x) x/sqrt(nCells), cellfun(@(x) std(x,[],2), fano_vis_cmlvCyc,'unif',false), 'unif',false) ;
fano_aud_cmlvCyc_ste = cellfun(@(x) x/sqrt(nCells), cellfun(@(x) std(x,[],2), fano_aud_cmlvCyc,'unif',false), 'unif',false) ;
fano_all_cmlvCyc_ste = cellfun(@(x) x/sqrt(nCells), cellfun(@(x) std(x,[],2), fano_all_cmlvCyc,'unif',false), 'unif',false) ;


base_win = [pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames];
resp_win = [trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames];

resp_vis_cmlvCyc_resp = cell(1,length(cycs));
resp_aud_cmlvCyc_resp= cell(1,length(cycs));
fano_vis_cmlvCyc_resp = cell(1,length(cycs));
fano_aud_cmlvCyc_resp = cell(1,length(cycs));
for icyc = 1:length(cycs)
    r = (trans_win(1):trans_win(end))+((icyc-1)*cycTime);
    b = (pre_win(1):pre_win(end))+((icyc-1)*cycTime);
   resp_vis_cmlvCyc_resp{icyc} = bsxfun(@minus, mean(resp_vis_cmlvCyc{icyc}(r,:),1), mean(resp_vis_cmlvCyc{icyc}(b,:),1)); 
   resp_aud_cmlvCyc_resp{icyc} = bsxfun(@minus, mean(resp_aud_cmlvCyc{icyc}(r,:),1), mean(resp_aud_cmlvCyc{icyc}(b,:),1));     
   fano_vis_cmlvCyc_resp{icyc} = bsxfun(@minus, mean(fano_vis_cmlvCyc{icyc}(r,:),1), mean(fano_vis_cmlvCyc{icyc}(b,:),1)); 
   fano_aud_cmlvCyc_resp{icyc} = bsxfun(@minus, mean(fano_aud_cmlvCyc{icyc}(r,:),1), mean(fano_aud_cmlvCyc{icyc}(b,:),1));     
end

pcmlvCyc = zeros(1,length(cycs));
pcmlvCycFano = zeros(1,length(cycs));
for icyc = 1:length(cycs)
[h, pcmlvCyc(icyc)] = ttest(resp_vis_cmlvCyc_resp{icyc}, resp_aud_cmlvCyc_resp{icyc},'alpha', 0.05/nCells);%, 'dim', 2, 'tail', 'left', 'alpha', 0.05/nCells);
[h, pcmlvCycFano(icyc)] = ttest(fano_vis_cmlvCyc_resp{icyc}, fano_aud_cmlvCyc_resp{icyc},'alpha', 0.05/nCells);%, 'dim', 2, 'tail', 'left', 'alpha', 0.05/nCells);
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

fanoAvsVFig_cmlvCycsErr = figure;
suptitle({'fano factor analysis';[titleStr '; alpha = ' num2str(0.05/nCells)]})
fanoAllTrialsFig_cmlvCycsErr = figure;
suptitle({'fano factor analysis';[titleStr '; alpha = ' num2str(0.05/nCells)]})
for icyc = 1:length(cycs)
figure(fanoAvsVFig_cmlvCycsErr);
subplot(c,c2,icyc)
shadedErrorBar(ttCycMs(1:size(fano_vis_cmlvCyc{icyc},1)), fano_vis_cmlvCyc_mean{icyc},fano_vis_cmlvCyc_ste{icyc},'g');
hold on
shadedErrorBar(ttCycMs(1:size(fano_aud_cmlvCyc{icyc},1)), fano_aud_cmlvCyc_mean{icyc},fano_aud_cmlvCyc_ste{icyc},'k');
hold on
vline((resp_win+((icyc-1)*cycTime)/(cycTime/cycTimeMs)), '--r')
vline((base_win+((icyc-1)*cycTime)/(cycTime/cycTimeMs)), '--k')
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
xlim([-10 cycTime*icyc]/(cycTime/cycTimeMs))
% ylim([-0.01 0.03]);
title([num2str(icyc) ' cycs;All cells; p = ' num2str(pcmlvCycFano(icyc))])

figure(fanoAllTrialsFig_cmlvCycsErr);
subplot(c,c2,icyc)
shadedErrorBar(ttCycMs(1:size(fano_all_cmlvCyc{icyc},1)), fano_all_cmlvCyc_mean{icyc},fano_vis_cmlvCyc_ste{icyc},'k');
hold on
vline((resp_win+((icyc-1)*cycTime)/(cycTime/cycTimeMs)), '--r')
vline((base_win+((icyc-1)*cycTime)/(cycTime/cycTimeMs)), '--k')
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
xlim([-10 cycTime*icyc]/(cycTime/cycTimeMs))
% ylim([-0.01 0.03]);
title([num2str(icyc) ' cycs;All cells'])
end

figure(respAvsVFig_cmlvCycs);
print([fnout 'press_align_TCbycmlvCyc' datasetStr '.pdf'], '-dpdf');
figure(respAvsVFig_cmlvCycsErr);
print([fnout 'press_align_TCbycmlvCycErr' datasetStr '.pdf'], '-dpdf')
figure(fanoAvsVFig_cmlvCycsErr);
print([fnout 'press_align_FanoCmlvCycErr_AvsV' datasetStr '.pdf'], '-dpdf');
figure(respAvsVFig_cmlvCycs);
print([fnout 'press_align_FanoCmlvCycErr_all' datasetStr '.pdf'], '-dpdf');

%% trials of a specific length
resp_vis_cyc = cell(1,length(cycs));
resp_aud_cyc = cell(1,length(cycs));
resp_all_cyc = cell(1,length(cycs));
nTrials_vis = zeros(1,length(cycs));
nTrials_aud = zeros(1,length(cycs));
endTrResp_vis = cell(1,length(cycs));
endTrResp_aud = cell(1,length(cycs));
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(respCellsInd).ind,cell_ind);
        if mouse(imouse).expt(iexp).info.cyc_time == 11
            for icyc = 1:length(cycs)
                if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp,2) >= icyc
            if size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{icyc},3) > 0 
                resp_vis_cyc{icyc} = cat(2,resp_vis_cyc{icyc},mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{icyc}(:,cell_ind,:),3));
                resp_aud_cyc{icyc} = cat(2,resp_aud_cyc{icyc},mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cycResp{icyc}(:,cell_ind,:),3));
                resp_all_cyc{icyc} = cat(2,resp_all_cyc{icyc},mean(cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{icyc}(:,cell_ind,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cycResp{icyc}(:,cell_ind,:)),3));
                nTrials_vis(icyc) = nTrials_vis(icyc)+size(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{icyc}(:,cell_ind,:),3);
                nTrials_aud(icyc) = nTrials_aud(icyc)+size(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cycResp{icyc}(:,cell_ind,:),3);
                endTrialRespIntervalFrames = (mouse(imouse).expt(iexp).info.cyc_time_ms - 100)*(mouse(imouse).expt(iexp).info.cyc_time/mouse(imouse).expt(iexp).info.cyc_time_ms)-1;
                endTrPreWin = mouse(imouse).expt(iexp).win(1).frames+((icyc-1)*mouse(imouse).expt(iexp).info.cyc_time);
                endTrResp_vis{icyc} = cat(2,endTrResp_vis{icyc},bsxfun(@minus,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{icyc}(end-endTrialRespIntervalFrames:end,cell_ind,:),3),1),mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cycResp{icyc}(endTrPreWin,cell_ind,:),3),1)));
                endTrResp_aud{icyc} = cat(2,endTrResp_aud{icyc},bsxfun(@minus,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cycResp{icyc}(end-endTrialRespIntervalFrames:end,cell_ind,:),3),1),mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cycResp{icyc}(endTrPreWin,cell_ind,:),3),1)));            

            end
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
    if ~isempty(resp_vis_cyc{icyc})
    r = (trans_win(1):trans_win(end))+((icyc-1)*cycTime);
    b = (pre_win(1):pre_win(end))+((icyc-1)*cycTime);
   resp_vis_cyc_resp{icyc} = bsxfun(@minus, mean(resp_vis_cyc{icyc}(r,:),1), mean(resp_vis_cyc{icyc}(b,:),1)); 
   resp_aud_cyc_resp{icyc} = bsxfun(@minus, mean(resp_aud_cyc{icyc}(r,:),1), mean(resp_aud_cyc{icyc}(b,:),1));     
    end
end

nCells = size(resp_vis_cyc_resp{6},2);

pCyc = zeros(1,length(cycs));
for icyc = 1:length(cycs)
    if ~isempty(resp_vis_cyc{icyc})
[h, pCyc(icyc)] = ttest(resp_vis_cyc_resp{icyc}, resp_aud_cyc_resp{icyc},'alpha', 0.05/nCells);%, 'dim', 2, 'tail', 'left', 'alpha', 0.05/nCells);
    end
end

resp_vis_cyc_mean = cellfun(@(x) nanmean(x,2),resp_vis_cyc,'unif',false);
resp_aud_cyc_mean = cellfun(@(x) nanmean(x,2),resp_aud_cyc,'unif',false);

ttCycMs = (-pre_event_frames:cycs(end)*cycTime)/(cycTime/cycTimeMs);
respAvsVFig_cycs = figure;
suptitle({[titleStr '; alpha = ' num2str(0.05/nCells)]})
for icyc = 1:length(cycs)
    if ~isempty(resp_vis_cyc{icyc})
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
end

resp_vis_cyc_ste = cellfun(@(x) x/sqrt(nCells), cellfun(@(x) nanstd(x,[],2), resp_vis_cyc,'unif',false), 'unif',false) ;
resp_aud_cyc_ste = cellfun(@(x) x/sqrt(nCells), cellfun(@(x) nanstd(x,[],2), resp_aud_cyc,'unif',false), 'unif',false) ;

respAvsVFig_cycsErr = figure;
suptitle({[titleStr '; alpha = ' num2str(0.05/nCells)]})
for icyc = 1:length(cycs)
    if ~isempty(resp_vis_cyc{icyc})
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
end

figure(respAvsVFig_cycs);
print([fnout 'press_align_TCbyCyc' datasetStr '.pdf'], '-dpdf');
figure(respAvsVFig_cycsErr);
print([fnout 'press_align_TCbyCycErr' datasetStr '.pdf'], '-dpdf')

%%
respByCycScatter = figure;
suptitle([titleStr '-resp from 100ms after last base stim to end'])
for icyc = 1:length(cycs)
    if ~isempty(endTrResp_vis{icyc})
        subplot(c,c2,icyc)
        scatter(endTrResp_vis{icyc},endTrResp_aud{icyc},50,'k.');
        hold on
        plot(-20:1:20,-20:1:20,'k--')
        ax = [min([endTrResp_vis{icyc} endTrResp_aud{icyc}]) max([endTrResp_vis{icyc} endTrResp_aud{icyc}])];
        xlim(ax)
        ylim(ax)
        axis square
        title([num2str(icyc) ' cycs;' num2str(nTrials_vis(icyc)) '/' num2str(nTrials_aud(icyc))])
    end
end
figure(respByCycScatter);
print([fnout 'press_align_EndTrialsRespbyCyc' datasetStr '.pdf'], '-dpdf');
end