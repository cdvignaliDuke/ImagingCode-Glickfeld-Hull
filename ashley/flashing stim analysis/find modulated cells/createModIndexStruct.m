function msModCells = createModIndexStruct(datasetStr)
%% load data, set data group
% % cellsInd = 14;
eval(['awFSAVdatasets' datasetStr])
titleStr = datasetStr;
if strcmp(titleStr, '')
    titleStr = 'V1';
else
    titleStr = titleStr(2:end);
end
rc = behavConstsAV;
if strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' datasetStr];
else
    dataGroup = [];
end
str = unique({expt.SubNum});
mouse_str = ['i' strjoin(str,'_i')];
% % if cellsOnly    
% % load(fullfile(rc.caOutputDir,dataGroup,'cells only', [mouse_str '_CaSummary' datasetStr '.mat']));
% % titleStr = [titleStr mouse(1).expt(1).cells(cellsInd).name];
% % fnout = fullfile(rc.caOutputDir,dataGroup,'cells only', [titleStr '_' mouse_str]);
% % else
load(fullfile(rc.caOutputDir,dataGroup, [mouse_str '_CaSummary' datasetStr '.mat']));
% % titleStr = [titleStr mouse(1).expt(1).cells(cellsInd).name];
fnout = fullfile(rc.caOutputDir,dataGroup, [titleStr '_' mouse_str]);
% % end
% % disp(fnout)
%% set params for figures
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
% set(0,'DefaultaxesFontSize', 16)
%% set analysis windows
earlyRespMaxMs = 1000;
lateRespMaxMs = 2800;

earlyWinMs = [500 1000];
lateWinMs = [2300 2800];

trialLengthCutoffMs = 1300;

trialDirectionCutoffDeg = 30;

nexp = 0;
for imouse = 1:size(mouse,2)
    nexp = nexp+size(mouse(imouse).expt,2);
end
n = ceil(sqrt(nexp+1));
if (n^2)-n > nexp+1
    n2 = n-1;
else
    n2 = n;
end

visual = 1;
auditory = 2;
catchAlign = 3;
targetAlign = 2;
entireTrAlign = 4;
startAlign = 1;

hits = 1;
misses = 2;
fas = 3;
crs = 4;

% for resp histograms for each cell
nbins = 20;
bin_edges = linspace(-0.2,0.2,nbins);

%% get experiment parameters
pre_win = mouse(1).expt(1).win(1).frames;
trans_win = mouse(1).expt(1).win(2).frames;
pre_event_frames = mouse(1).expt(1).pre_event_frames;
post_event_frames = mouse(1).expt(1).post_event_frames;

tt = -pre_event_frames:post_event_frames-1;
% ttMs = tt/(cycTime/cycTimeMs);
% baseStimFrames = -(floor(pre_event_frames/cycTime)*cycTime):cycTime:0;

%% find responses of each neuron for each modulation analysis
iplot = 1;
AVcdfFig = figure;
suptitle('response to target stim compared to proximal baseline vis stim')
val_inv_cdfFig = figure;
suptitle('response to target stim compared to proximal baseline vis stim')
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)  
        %experiment specifics, so that both 100ms and 200ms stim on trials
        %can be used
        cycTime = mouse(imouse).expt(iexp).info(1).cyc_time;
        cycTimeMs = mouse(imouse).expt(iexp).info(1).cyc_time_ms;
        earlyCyc = ceil(earlyRespMaxMs/cycTimeMs);
        lateCyc = ceil(lateRespMaxMs/cycTimeMs);
        trialLengthCutoffCyc = ceil(trialLengthCutoffMs/cycTimeMs);
        earlyWinFrames = earlyWinMs*(cycTime/cycTimeMs)+pre_event_frames;
        lateWinFrames = lateWinMs*(cycTime/cycTimeMs)+pre_event_frames;
        for iav = 1:2 %for both visual and auditory target enhancement 
            for iout = 1:4 %for hits, misses, invalid hits, invalid misses
                if iout > 2 % only do iout if the experiment included catch trials and iav is on visual trials (1)
                    if mouse(imouse).expt(iexp).info.isCatch == 0
                        continue
                    end
                    if iav == 2
                        continue
                    end
                end
                
                tc_cyc = mouse(imouse).expt(iexp).align(entireTrAlign).av(iav).outcome(iout).cycResp; %aligned to trial start for each cycle length
                cycs = num2cell(1:length(tc_cyc));
                target = mouse(imouse).expt(iexp).align(entireTrAlign).av(iav).outcome(iout).direction;
                target_small = cellfun(@(x) find(x < trialDirectionCutoffDeg),target,'unif',false);
                target_large = cellfun(@(x) find(x > trialDirectionCutoffDeg),target,'unif',false);
                
                
                shortCycs = intersect(find(cell2mat(cellfun(@(x) ~isempty(x),tc_cyc,'unif',false))),1:trialLengthCutoffCyc); % sort cycles into short and long types
                longCycs = intersect(find(cell2mat(cellfun(@(x) ~isempty(x),tc_cyc,'unif',false))),trialLengthCutoffCyc+1:length(tc_cyc));

                pre_base1 = cellfun(@(x) squeeze(mean(x(pre_win,:,:),1)),tc_cyc(2:end),'unif',false); % find response to last baseline vis stim
                trans_base1 = cellfun(@(x) squeeze(mean(x(trans_win,:,:),1)),tc_cyc(2:end),'unif',false); 
                base1_cyc = [{[]}, cellfun(@(x,y) bsxfun(@minus, x,y), trans_base1, pre_base1,'unif',false)];
                
                
                pre_base = cellfun(@(x,y) squeeze(mean(x(pre_win+(cycTime*(y-1)),:,:),1)),tc_cyc(2:end),cycs(2:end),'unif',false); % find response to first baseline vis stim
                trans_base = cellfun(@(x,y) squeeze(mean(x(trans_win+(cycTime*(y-1)),:,:),1)),tc_cyc(2:end),cycs(2:end),'unif',false); 
                base_cyc = [{[]}, cellfun(@(x,y) bsxfun(@minus, x,y), trans_base, pre_base,'unif',false)];

                pre_tar = cellfun(@(x,y) squeeze(mean(x(pre_win+(cycTime*(y)),:,:),1)),tc_cyc(2:end),cycs(2:end),'unif',false);% find response to target for each trial
                trans_tar = cellfun(@(x,y) squeeze(mean(x(trans_win+(cycTime*(y)),:,:),1)),tc_cyc(2:end),cycs(2:end),'unif',false);
                tar_cyc = [{[]}, cellfun(@(x,y) bsxfun(@minus, x,y), trans_tar, pre_tar,'unif',false)];
                
                
                ind = find(cell2mat(cellfun(@(x) length(size(x)),tc_cyc,'unif',false)) == 2);% need to do this to correct for reshape that happens when there is only one trial for that cycle type
                base_cyc(ind) = cellfun(@transpose,base_cyc(ind),'unif',false);
                base1_cyc(ind) = cellfun(@transpose,base1_cyc(ind),'unif',false);
                tar_cyc(ind) = cellfun(@transpose,tar_cyc(ind),'unif',false);
                
                if sum(cell2mat(cellfun(@(x) size(x,3),tc_cyc(shortCycs),'unif',false))) > 1 %combine response to short and long trials, only if there is more than one total trial
%                     if any (cell2mat(cellfun(@(x) length(size(x)),tc_cyc(shortCycs),'unif',false)) < 3) 
%                         ind = find(cell2mat(cellfun(@(x) length(size(x)),tc_cyc(shortCycs),'unif',false)) == 2);
%                         base_cyc(shortCycs(ind)) = cellfun(@transpose,base_cyc(shortCycs(ind)),'unif',false);
%                         tar_cyc(shortCycs(ind)) = cellfun(@transpose,tar_cyc(shortCycs(ind)),'unif',false);
%                     end
                        base_cyc_short =  cat(2,cell2mat(base_cyc(shortCycs)));
                        tar_cyc_short = cat(2,cell2mat(tar_cyc(shortCycs)));
                else                        
                    base_cyc_short = [];
                    tar_cyc_short = [];
                end
                if sum(cell2mat(cellfun(@(x) size(x,3),tc_cyc(longCycs),'unif',false))) > 1  
%                     if any (cell2mat(cellfun(@(x) length(size(x)),tc_cyc(longCycs),'unif',false)) < 3)
%                         base_cyc(longCycs(ind)) = cellfun(@transpose,base_cyc(longCycs(ind)),'unif',false);
%                         tar_cyc(longCycs(ind)) = cellfun(@transpose,tar_cyc(longCycs(ind)),'unif',false);
%                     end
                        base_cyc_long = cat(2,cell2mat(base_cyc(longCycs)));
                        tar_cyc_long = cat(2,cell2mat(tar_cyc(longCycs)));
                else                        
                    base_cyc_long = [];
                    tar_cyc_long = [];
                end
                
                if sum(cell2mat(cellfun(@(x) size(x,3),tc_cyc,'unif',false))) > 1  
%                     if any (cell2mat(cellfun(@(x) length(size(x)),tc_cyc,'unif',false)) < 3)
%                         ind = find(cell2mat(cellfun(@(x) length(size(x)),tc_cyc,'unif',false)) == 2);
%                         base_cyc(ind) = cellfun(@transpose,base_cyc(ind),'unif',false);
%                         tar_cyc(ind) = cellfun(@transpose,tar_cyc(ind),'unif',false);
%                     end
                        ind = find(cell2mat(cellfun(@(x) ~isempty(x),base_cyc,'unif',false)));
                        base_cyc_all = cat(2,cell2mat(base_cyc));
                        base1_cyc_all = cat(2,cell2mat(base1_cyc));
                        tar_cyc_all = cat(2,cell2mat(tar_cyc));
                        target_all = chop(cat(2,cell2mat(target(ind))),2);
                else                        
                    base_cyc_all = [];
                    base1_cyc_all = [];
                    tar_cyc_all = [];
                    target_all = [];
                end
                
                if iav == 1 
                    if iout < 3 % for hits and misses, find response to visual and auditory trials, early and late in the trial
                    tc_early_v = mouse(imouse).expt(iexp).align(startAlign).av(visual).outcome(iout).cmlvCycResp{earlyCyc};
                    tc_late_v = mouse(imouse).expt(iexp).align(startAlign).av(visual).outcome(iout).cmlvCycResp{lateCyc};
                    tc_early_a = mouse(imouse).expt(iexp).align(startAlign).av(auditory).outcome(iout).cmlvCycResp{earlyCyc};
                    tc_late_a = mouse(imouse).expt(iexp).align(startAlign).av(auditory).outcome(iout).cmlvCycResp{lateCyc};

                    if size(tc_early_v,3) > 1 & size(tc_early_a,3) > 1 
                        tc_early_v_resp = squeeze(mean(tc_early_v(earlyWinFrames,:,:),1));
                        tc_early_a_resp = squeeze(mean(tc_early_a(earlyWinFrames,:,:),1));
                    else                        
                        tc_early_v_resp = [];
                        tc_early_a_resp = [];
                    end
                    if size(tc_late_v,3) > 1 & size(tc_late_a,3) > 1 
                        tc_late_v_resp = squeeze(mean(tc_late_v(lateWinFrames,:,:),1));
                        tc_late_a_resp = squeeze(mean(tc_late_a(lateWinFrames,:,:),1));  
                    else
                        tc_late_v_resp = [];
                        tc_late_a_resp = [];              
                    end
                    % find modulation index (V-A)/(V+A) for early and late
                    % responses
                    early_resp_avmi = bsxfun(@rdivide, bsxfun(@minus, mean(tc_early_v_resp,2) , mean(tc_early_a_resp,2)), bsxfun(@plus, mean(tc_early_v_resp,2) , mean(tc_early_a_resp,2)));
                    late_resp_avmi = bsxfun(@rdivide, bsxfun(@minus, mean(tc_late_v_resp,2) , mean(tc_late_a_resp,2)), bsxfun(@plus, mean(tc_late_v_resp,2) , mean(tc_late_a_resp,2)));
                    end
                    if mouse(imouse).expt(iexp).info.isCatch
                        if iout == 1
                            tc_target_val = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp;    
                            tc_target_inv = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp;
                            tar_resp_val = cellfun(@minus, cellfun(@(x) mean(mean(x(pre_win,:,:),3),1),tc_target_val,'unif',false), cellfun(@(x) mean(mean(x(trans_win,:,:),3),1),tc_target_val,'unif',false),'unif',false);
                            tar_resp_inv = cellfun(@minus, cellfun(@(x) mean(mean(x(pre_win,:,:),3),1),tc_target_inv,'unif',false), cellfun(@(x) mean(mean(x(trans_win,:,:),3),1),tc_target_inv,'unif',false),'unif',false);
                            % find modulation index (val-inv)/(val+inv)
                            topmi = cellfun(@(x,y) bsxfun(@minus,x,y),tar_resp_val,tar_resp_inv,'unif',false);
                            botmi = cellfun(@(x,y) bsxfun(@plus,x,y),tar_resp_val,tar_resp_inv,'unif',false);
                            tar_resp_mi = cellfun(@(x,y) bsxfun(@rdivide,x,y),topmi,botmi,'unif',false);
                        end
                    end
                end
                ncells = size(tc_early_v,2);

                for imi = 1:2 

                    if imi == 1 %compare trial types by area under ROC curve, for each trial outcome type
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).name = 'roc auc';   
                    % area under ROC curve where noise is response to
                    % baseline, signal is response to target, short and
                    % long trials
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).name = 'target enhancement over last base stim';
                        if ~isempty(base_cyc_short) & ~isempty(base_cyc_short)
                        auc_short = nan(ncells,1);
                        rs_short = nan(ncells,1);
                        for icell = 1:ncells
                            auc_short(icell) = roc_gh(base_cyc_short(icell,:),tar_cyc_short(icell,:));
                            [rs_tar_p rs_short(icell)] = ranksum(base_cyc_short(icell,:),tar_cyc_short(icell,:));
                        end
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(1).name = 'short trials';
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(1).value = auc_short;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(1).rst = find(rs_short);
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(1).ind = find(auc_short > 0.5);
                        end                        
                       
                        if ~isempty(base_cyc_long) & ~isempty(tar_cyc_long)
                        auc_long = nan(1,ncells);
                        rs_long = nan(1,ncells);
                        for icell = 1:ncells
                            auc_long(icell) = roc_gh(base_cyc_long(icell,:),tar_cyc_long(icell,:));
                            [rs_tar_p rs_long(icell)] = ranksum(base_cyc_long(icell,:),tar_cyc_long(icell,:));
                        end
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(2).name = 'long trials';
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(2).value = auc_long;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(2).rst = find(rs_long);
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(2).ind = find(auc_long > 0.5);
                        end
                        
                        if ~isempty(base_cyc_all) & ~isempty(tar_cyc_all)
                        auc_all = nan(1,ncells);
                        rs_all = nan(1,ncells);
                        rs_p = nan(1,ncells);
                        resp_ratio = nan(1,ncells);
%                         rh_base = nan(ncells,nbins);
%                         rh_tar = nan(ncells,nbins);
                        dirs = unique(target_all);
                        auc_tar_all = cell(1,length(dirs));
                        rs_tar_all = cell(1,length(dirs));
%                         resp_base_all = cell(1,length(dirs));
%                         resp_tar_all = cell(1,length(dirs));
                        resp_ratio_dirs = cell(1,length(dirs));
                        for icell = 1:ncells
                            resp_ratio(icell) = mean(tar_cyc_all(icell,:))/mean(base_cyc_all(icell,:));
                            auc_all(icell) = roc_gh(base_cyc_all(icell,:),tar_cyc_all(icell,:));
                            [rs_p(icell) rs_all(icell)] = ranksum(base_cyc_all(icell,:),tar_cyc_all(icell,:));
%                             rh_base(icell,:) = histc(base_cyc_all(icell,:),bin_edges);
%                             rh_tar(icell,:) = histc(tar_cyc_all(icell,:),bin_edges);
                            for idir = 1:length(dirs)
                                ind = find(target_all == dirs(idir));
                                auc_tar_all{idir} = cat(1,auc_tar_all{idir},roc_gh(base_cyc_all(icell,ind),tar_cyc_all(icell,ind)));
                                resp_ratio_dirs{idir} = cat(1,resp_ratio_dirs{idir},mean(tar_cyc_all(icell,ind))/mean(base_cyc_all(icell,ind)));
                                
                                [rs_tar_p rs_temp] = ranksum(base_cyc_all(icell,ind),tar_cyc_all(icell,ind));
                                rs_tar_all{idir} = cat(1,rs_tar_all{idir},rs_temp);
%                                 resp_base_all{idir} = cat(1,resp_base_all{idir},histc(base_cyc_all(icell,ind),bin_edges));
%                                 resp_tar_all{idir} = cat(1,resp_tar_all{idir},histc(tar_cyc_all(icell,ind),bin_edges));
                                
                            end
                            
                        end
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).name = 'all trials';
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).value = auc_all;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).rst_p = rs_p;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).rst = rs_all;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).ratio = resp_ratio;
%                         msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).resp_base = rh_base;
%                         msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).resp_tar = rh_tar;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).ind = find(auc_all > 0.5);
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).dirs = dirs;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).value_dirs = auc_tar_all;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).rst_dirs_logical = rs_tar_all;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).rst_dirs = cellfun(@(x) find(x),rs_tar_all,'unif',false);
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).ratio_dirs = resp_ratio_dirs;
%                         msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).resp_base_dirs = resp_base_all;
%                         msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).resp_tar_dirs = resp_tar_all;
                        
                        
                        %compare target to first base stim resp
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).name = 'target enhancement over first base stim';
                        auc_all = nan(1,ncells);
                        rs_all = nan(1,ncells);
                        auc_tar_all = cell(1,length(dirs));
                        rs_tar_all = cell(1,length(dirs));
                        for icell = 1:ncells
                            auc_all(icell) = roc_gh(base1_cyc_all(icell,:),tar_cyc_all(icell,:));
                            [rs_p rs_all(icell)] = ranksum(base1_cyc_all(icell,:),tar_cyc_all(icell,:));
                            
                            for idir = 1:length(dirs)
                                ind = find(target_all == dirs(idir));
                                auc_tar_all{idir} = cat(1,auc_tar_all{idir},roc_gh(base1_cyc_all(icell,ind),tar_cyc_all(icell,ind)));
                                [rs_tar_p rs_temp] = ranksum(base1_cyc_all(icell,ind),tar_cyc_all(icell,ind));
                                rs_tar_all{idir} = cat(1,rs_tar_all{idir},rs_temp);
                                
                            end
                            
                        end
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).name = 'all trials';
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).value = auc_all;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).ustat = rs_all;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).ind = find(auc_all > 0.5);
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).dirs = dirs;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).value_dirs = auc_tar_all;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).rst_dirs = cellfun(@(x) find(x),rs_tar_all,'unif',false);
                        
                        else
                            
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).name = 'all trials';
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).value = [];
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).rst_p = [];
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).rst = [];
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).ind = [];
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).dirs = [];
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).value_dirs = [];
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(3).rst_dirs = [];
                            
                        
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).name = 'all trials';
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).value = [];
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).ustat = [];
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).ind =[];
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).dirs = [];
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).value_dirs = [];
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(3).auc(3).rst_dirs = [];
                        end

                        if iav == 2 & iout == 1
                            figure(AVcdfFig);
                            subplot(n,n2,iplot)
                            AVcdf(msModCells(imouse).expt(iexp).av(1).outcome(iout).mi(imi).comp(1).auc(3).value,msModCells(imouse).expt(iexp).av(2).outcome(iout).mi(imi).comp(1).auc(3).value,[mouse(imouse).expt(iexp).mouse_name '-' mouse(imouse).expt(iexp).date])
                            xlabel('auc')
                        end
                        
                        if iav == 1
                            
                            %ROC analysis where auditory trial is noise,
                            %visual trial is signal
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(2).name = 'av_resp';
                        if iout < 3
                        if ~isempty(tc_late_a_resp) & ~isempty(tc_late_v_resp)
                        auc_late = nan(1,ncells);
                        for icell = 1:ncells
                            auc_late(icell) = roc_gh(tc_late_a_resp(icell,:),tc_late_v_resp(icell,:));
                        end
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(2).auc(2).name = 'late phase';
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(2).auc(2).value = auc_late;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(2).auc(2).ind = find(auc_late > 0.5);
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(3).comp(1).sub(2).name = 'late phase';
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(3).comp(1).sub(2).value = late_resp_avmi;
                        end
                        if  ~isempty(tc_early_a_resp) & ~isempty(tc_early_v_resp)
                        auc_early = nan(1,ncells);
                        for icell = 1:ncells
                            auc_early(icell) = roc_gh(tc_early_a_resp(icell,:),tc_early_v_resp(icell,:));
                        end
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(2).auc(1).name = 'early phase';
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(2).auc(1).value = auc_early;
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(2).auc(1).ind = find(auc_early > 0.5);
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(3).comp(1).sub(1).name = 'early phase';
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(3).comp(1).sub(1).value = early_resp_avmi;
                        end

                        % normalized modulation indices, across visual and
                        % auditory trial types, and valid and invalid
                        % visual trials
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(3).name = 'norm sub';    
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(3).comp(1).name = '(V-A)/(V+A)'; 
                        if iout == 1 & mouse(imouse).expt(iexp).info.isCatch
                            msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(3).comp(2).name = '(val-inv)/(val+inv)';
                            msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(3).comp(2).value = tar_resp_mi;
                        end
                        end
                        
                        end
                    elseif imi == 2 % compare trial types by ttest
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).name = 'ttest';

                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).name = 'target enhancement';
                        ttest_short = find(ttest(base_cyc_short, tar_cyc_short, 'tail', 'left'));
                        ttest_long = find(ttest(base_cyc_long, tar_cyc_long, 'tail', 'left'));

                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(1).name = 'short trials';
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(1).ind = ttest_short;

                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(2).name = 'long trials';
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(1).auc(2).ind = ttest_long;

                        if iav == 1
                        msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(2).name = 'av_resp';

%                         ttest_early = find(ttest(tc_early_a_resp,tc_early_v_resp,'dim',1, 'tail', 'left'));
%                         ttest_late = find(ttest(tc_late_a_resp,tc_late_v_resp,'dim',1, 'tail', 'left')); % not sure this is the right variable/calc??? need same number of trials in each?
% 
%                         msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(2).auc(1).name = 'early phase';
%                         msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(2).auc(1).ind = ttest_early;
%                         msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(2).auc(2).name = 'late phase';
%                         msModCells(imouse).expt(iexp).av(iav).outcome(iout).mi(imi).comp(2).auc(2).ind = ttest_late;

                        end        
                    end
                end
            end
        end
       iplot = iplot+1; 
    end
end



%% save analysis
    save([fnout '_modCells'],'msModCells');
end