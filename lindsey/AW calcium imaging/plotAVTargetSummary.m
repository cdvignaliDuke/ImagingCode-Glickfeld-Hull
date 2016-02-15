function plotAVTargetSummary(datasetStr)
%takes mouse structure and plots transient responses to target and base
%stimuli
close all
av = behavParamsAV;
eval(['awFSAVdatasets' datasetStr])
titleStr = datasetStr;
if strcmp(titleStr, '')
    titleStr = 'V1';
else
    titleStr = titleStr(2:end);
end
rc = behavConstsAV;
str = unique({expt.SubNum});
values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
mouse_str = [];
for imouse = 1:size(str,2)
    mouse_str = [mouse_str 'i' str{imouse} '_'];  
end
load(fullfile(rc.caOutputDir, [mouse_str 'CaSummary' datasetStr '.mat']));
pre_win = mouse(1).expt(1).win(1).frames;
trans_win = mouse(1).expt(1).win(2).frames;
pre_event_frames = mouse(1).expt(1).pre_event_frames;
post_event_frames = mouse(1).expt(1).post_event_frames;
fnout = fullfile(rc.caOutputDir, [date '_' mouse_str]);

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
figure;
i = 1;
tt = -pre_event_frames:post_event_frames-1;
resp = [];
resp_aud = [];
resp_FA = [];
resp_CR = [];
ialign = 2;
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        subplot(3,4,i)
        i = i+1;
        good_ind = mouse(imouse).expt(iexp).cells(1).ind;
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(3).resp(:,good_ind,:),2),3), 'c');
        hold on
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).resp(:,good_ind,:),2),3), 'k');
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).resp(:,good_ind,:),2),3), 'g');
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(4).resp(:,good_ind,:),2),3), 'b');
        vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
        vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
        xlim([-10 30])
        resp = cat(2, resp, mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).resp(:,good_ind,:),3));
        resp_aud = cat(2, resp_aud, mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).resp(:,good_ind,:),3));
        resp_FA = cat(2, resp_FA, mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(3).resp(:,good_ind,:),3));
        resp_CR = cat(2, resp_CR, mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(4).resp(:,good_ind,:),3));
        title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(good_ind)) ' cells']})
    end
end
subplot(3,4,i)
plot(tt, nanmean(resp,2), 'k')
hold on
plot(tt, nanmean(resp_aud,2), 'g')
plot(tt, nanmean(resp_FA,2), 'c')
plot(tt, nanmean(resp_CR,2), 'b')
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
xlim([-10 30])
title(['All cells; n = ' num2str(size(resp,2))])
suptitle({titleStr, 'Hits: Black; FAs: Cyan; CR: Blue; Auditory: Green'})
print([fnout 'release_align_TCs' datasetStr '.pdf'], '-dpdf')



figure;
i = 1;
tt = -pre_event_frames:post_event_frames-1;
resp = [];
resp_aud = [];
ialign = 1;
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        subplot(3,4,i)
        i = i+1;
        good_ind = mouse(imouse).expt(iexp).cells(1).ind;
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).resp(:,good_ind,:),2),3), 'k');
        hold on
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).resp(:,good_ind,:),2),3), 'g');
        vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
        vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
        xlim([-10 30])
        resp = cat(2,resp, mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).resp(:,good_ind,:),3));
        resp_aud = cat(2,resp_aud, mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).resp(:,good_ind,:),3));
        title({mouse(imouse).expt(iexp).date, ['n = ' num2str(length(good_ind)) ' cells']})
    end
end
subplot(3,4,i)
plot(tt,nanmean(resp,2), 'k')
hold on
plot(tt,nanmean(resp_aud,2), 'g')
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
xlim([-10 30])
title(['All cells; n = ' num2str(size(resp,2))])
suptitle({titleStr, 'Hits: Black; Auditory: Green'})
print([fnout 'press_align_TCs' datasetStr '.pdf'], '-dpdf')


for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        for ialign = 2
            for iav = 1:2
                for io = [1 2 5]
                    for iDir = 1:size(mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).stimResp,2)
                        tempresp = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).stimResp{iDir};
                        trans_resp = mean(tempresp(trans_win,:,:),1)-mean(tempresp(pre_win,:,:),1);
                        mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).trans_resp(:,iDir) =  squeeze(mean(trans_resp,3));
                    end
                end
                for io = [3 4]
                    tempresp = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).resp;
                    trans_resp = mean(tempresp(trans_win,:,:),1)-mean(tempresp(pre_win,:,:),1);
                    mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).trans_resp(1,:) =  squeeze(mean(trans_resp,3));
                end
            end
        end
    end
end
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        for ialign = 1
            for iav = 1:2
                for io = [1 2 5]
                	tempresp = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).resp;
                    trans_resp = mean(tempresp(trans_win,:,:),1)-mean(tempresp(pre_win,:,:),1);
                    mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).trans_resp(1,:) =  squeeze(mean(trans_resp,3));
                end
            end
        end
    end
end

base_dirs_all = [];
base_dirs_sub = [];
bd = 0;
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        if length(mouse(imouse).expt(iexp).cells(1).ind) > 10
            base_dirs_all = [base_dirs_all chop(mouse(imouse).expt(iexp).visTargets,2)];
            bd = bd+1;
        else
            base_dirs_sub = [base_dirs_sub chop(mouse(imouse).expt(iexp).visTargets,2)];
        end
    end
end
if bd./size(expt,2) < 0.25
    bd_u = unique(base_dirs_sub);
    counts = hist(base_dirs_sub, bd_u)./(size(expt,2)-bd);
    base_dirs_all = [base_dirs_all bd_u(find(counts>=0.5))];
end
if bd == 0;
    for imouse = 1:size(mouse,2)
        for iexp = 1:size(mouse(imouse).expt,2)
        	base_dirs_all = [base_dirs_all chop(mouse(imouse).expt(iexp).visTargets,2)];
        end
    end
    base_dirs = unique(base_dirs_all);
else
    base_dirs_each = unique(base_dirs_all);
    if size(expt,2) > 4
        base_count = zeros(1, length(base_dirs_each));
        for idir = 1:length(base_dirs_each)
            base_count(:,idir) = length(find(base_dirs_all == base_dirs_each(idir)));
        end
        base_dirs = base_dirs_each(find((base_count./max(base_count,[],2))>=0.5));
    else
        base_dirs = base_dirs_each;
    end
end

trans_target_resp_all = [];
trans_target_resp_notori = [];
trans_base_resp_all = cell(1,2);
trans_base_resp_notori = cell(1,2);
trans_target_resp = cell(1,4);
trans_target_resp_SM = cell(2,4);
trans_target_resp_all_SM = cell(1,2);
trans_target_resp_notori_SM = cell(1,2);
trans_base_resp_all_FC = cell(1,2);
trans_base_resp_notori_FC = cell(1,2);
trans_base_resp = cell(2,4);
trans_base_resp_FC = cell(2,4);
ncells = cell(size(mouse,2),4);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        good_ind = mouse(imouse).expt(iexp).cells(1).ind;
        if sum(ismember(chop(mouse(imouse).expt(iexp).visTargets,2), base_dirs),2)==length(mouse(imouse).expt(iexp).visTargets)
            if sum(ismember(base_dirs, chop(mouse(imouse).expt(iexp).visTargets,2)),2)==length(base_dirs)
                dir_ind = 1:length(base_dirs);
                isdir = ones(1,length(base_dirs));
            else
              [isdir dir_ind] = ismember(base_dirs, chop(mouse(imouse).expt(iexp).visTargets,2));
              dir_ind(dir_ind==0)=[];  
            end
        else
            [isdir dir_ind] = ismember(base_dirs, chop(mouse(imouse).expt(iexp).visTargets,2));
            dir_ind(dir_ind==0)=[];
        end
        all_ori_ind = [];
        for iOri = 1:4
            if iexp == 1
                ncells{imouse, iOri} = 0;
            end
            ori_ind = intersect(good_ind, mouse(imouse).expt(iexp).cells(iOri+1).ind);
            all_ori_ind = [all_ori_ind; ori_ind];
            ncells{imouse, iOri} = ncells{imouse, iOri} + length(ori_ind);
            if length(ori_ind)>0
                tempR = nan(length(ori_ind), length(base_dirs));
                tempR(:,find(isdir)) = mouse(imouse).expt(iexp).align(2).av(1).outcome(1).trans_resp(ori_ind,dir_ind);
                trans_target_resp{iOri} = [trans_target_resp{iOri}; tempR];
                tempR = nan(length(ori_ind), length(base_dirs));
                tempR(:,find(isdir)) = mouse(imouse).expt(iexp).align(2).av(1).outcome(5).trans_resp(ori_ind,dir_ind);
                trans_target_resp_SM{1,iOri} = [trans_target_resp_SM{1,iOri}; tempR];
                tempR = nan(length(ori_ind), length(base_dirs));
                tempR(:,find(isdir)) = mouse(imouse).expt(iexp).align(2).av(1).outcome(2).trans_resp(ori_ind,dir_ind);
                trans_target_resp_SM{2,iOri} = [trans_target_resp_SM{2,iOri}; tempR];
                for av = 1:2
                    trans_base_resp{av,iOri} = [trans_base_resp{av,iOri}; mouse(imouse).expt(iexp).align(1).av(av).outcome(1).trans_resp(1,ori_ind)'];
                end
                trans_base_resp_FC{1,iOri} = [trans_base_resp_FC{1,iOri}; mouse(imouse).expt(iexp).align(2).av(1).outcome(3).trans_resp(1,ori_ind)'];
                trans_base_resp_FC{2,iOri} = [trans_base_resp_FC{2,iOri}; mouse(imouse).expt(iexp).align(2).av(1).outcome(4).trans_resp(1,ori_ind)'];
            elseif iexp == 1 & imouse == 1
                trans_target_resp{iOri} = nan(1, length(base_dirs));
                for i = 1:2
                    trans_target_resp_SM{i,iOri} = nan(1, length(base_dirs));
                    trans_base_resp{i,iOri} = nan(1, 1);
                    trans_base_resp_FC{i,iOri} = nan(1, 1);
                end
            end
        end
        if length(all_ori_ind)>0
            tempR_all = nan(length(all_ori_ind), length(base_dirs));
            tempR_all(:,find(isdir)) = mouse(imouse).expt(iexp).align(2).av(1).outcome(1).trans_resp(all_ori_ind,dir_ind);
            trans_target_resp_all = [trans_target_resp_all; tempR_all];
            tempR_all = nan(length(all_ori_ind), length(base_dirs));
            tempR_all(:,find(isdir)) = mouse(imouse).expt(iexp).align(2).av(1).outcome(5).trans_resp(all_ori_ind,dir_ind);
            trans_target_resp_all_SM{1} = [trans_target_resp_all_SM{1}; tempR_all];
            tempR_all = nan(length(all_ori_ind), length(base_dirs));
            tempR_all(:,find(isdir)) = mouse(imouse).expt(iexp).align(2).av(1).outcome(2).trans_resp(all_ori_ind,dir_ind);
            trans_target_resp_all_SM{2} = [trans_target_resp_all_SM{2}; tempR_all];
            trans_base_resp_all_FC{1} = [trans_base_resp_all_FC{1}; mouse(imouse).expt(iexp).align(2).av(1).outcome(3).trans_resp(:,all_ori_ind)']; 
            trans_base_resp_all_FC{2} = [trans_base_resp_all_FC{2}; mouse(imouse).expt(iexp).align(2).av(1).outcome(4).trans_resp(:,all_ori_ind)'];
            for iav = 1:2
                trans_base_resp_all{iav} = [trans_base_resp_all{iav}; mouse(imouse).expt(iexp).align(1).av(iav).outcome(1).trans_resp(1,all_ori_ind)'];
            end
        elseif iexp == 1 & imouse == 1
            trans_target_resp_all = nan(1, length(base_dirs));
            for i = 1:2
                trans_target_resp_all_SM{i} = nan(1, length(base_dirs));
                trans_base_resp_all{i} = nan(1, 1);
                trans_base_resp_all_FC{i} = nan(1, 1);
            end
        end
        notori_ind = setdiff(good_ind, all_ori_ind);
        if length(notori_ind)>0
            tempR_notori = nan(length(notori_ind), length(base_dirs));
            tempR_notori(:,find(isdir)) = mouse(imouse).expt(iexp).align(2).av(1).outcome(1).trans_resp(notori_ind,dir_ind);
            trans_target_resp_notori = [trans_target_resp_notori; tempR_notori];
            tempR_notori = nan(length(notori_ind), length(base_dirs));
            tempR_notori(:,find(isdir)) = mouse(imouse).expt(iexp).align(2).av(1).outcome(5).trans_resp(notori_ind,dir_ind);
            trans_target_resp_notori_SM{1} = [trans_target_resp_notori_SM{1}; tempR_notori];
            tempR_notori = nan(length(notori_ind), length(base_dirs));
            tempR_notori(:,find(isdir)) = mouse(imouse).expt(iexp).align(2).av(1).outcome(2).trans_resp(notori_ind,dir_ind);
            trans_target_resp_notori_SM{2} = [trans_target_resp_notori_SM{2}; tempR_notori];
            trans_base_resp_notori_FC{1} = [trans_base_resp_notori_FC{1}; mouse(imouse).expt(iexp).align(2).av(1).outcome(3).trans_resp(:,notori_ind)']; 
            trans_base_resp_notori_FC{2} = [trans_base_resp_notori_FC{2}; mouse(imouse).expt(iexp).align(2).av(1).outcome(4).trans_resp(:,notori_ind)'];
            for iav = 1:2
                trans_base_resp_notori{iav} = [trans_base_resp_notori{iav}; mouse(imouse).expt(iexp).align(1).av(iav).outcome(1).trans_resp(1,notori_ind)'];
            end
        elseif iexp == 1 & imouse == 1
            trans_target_resp_notori = nan(1, length(base_dirs));
            for i = 1:2
                trans_target_resp_notori_SM{i} = nan(1, length(base_dirs));
                trans_base_resp_notori{i} = nan(1, 1);
                trans_base_resp_notori_FC{i} = nan(1, 1);
            end
        end
    end
end

figure;
for iOri = 1:4
    subplot(3,2,iOri)
    errorbar(base_dirs, nanmean(trans_target_resp{iOri},1), nanstd(trans_target_resp{iOri},[],1)./sqrt(sum(~isnan(trans_target_resp{iOri}),1)), '-ob');
    hold on
    errorbar(0, nanmean(trans_base_resp{1,iOri},1), nanstd(trans_base_resp{1,iOri},[],1)./sqrt(sum(~isnan(trans_base_resp{1,iOri}),1)), 'ok');
    errorbar(0, nanmean(trans_base_resp{2,iOri},1), nanstd(trans_base_resp{2,iOri},[],1)./sqrt(sum(~isnan(trans_base_resp{2,iOri}),1)), 'og');
    title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp{iOri},2)),1))])
    xlim([-10 100])
    hline(0, '--k')
end
subplot(3,2,5)
errorbar(base_dirs, nanmean(trans_target_resp_notori,1), nanstd(trans_target_resp_notori,[],1)./sqrt(sum(~isnan(trans_target_resp_notori),1)), '-ob');
hold on
errorbar(0, nanmean(trans_base_resp_notori{1},1), nanstd(trans_base_resp_notori{1},[],1)./sqrt(sum(~isnan(trans_base_resp_notori{1}),1)), 'ok');
errorbar(0, nanmean(trans_base_resp_notori{2},1), nanstd(trans_base_resp_notori{2},[],1)./sqrt(sum(~isnan(trans_base_resp_notori{2}),1)), 'og');
xlim([-10 100])
hline(0, '--k')
title(['Untuned- n = ' num2str(sum(~isnan(nanmean(trans_target_resp_notori,2)),1))])
subplot(3,2,6)
errorbar(base_dirs, nanmean(trans_target_resp_all,1), nanstd(trans_target_resp_all,[],1)./sqrt(sum(~isnan(trans_target_resp_all),1)), '-ob');
hold on
errorbar(0, nanmean(trans_base_resp_all{1},1), nanstd(trans_base_resp_all{1},[],1)./sqrt(sum(~isnan(trans_base_resp_all{1}),1)), 'ok');
errorbar(0, nanmean(trans_base_resp_all{2},1), nanstd(trans_base_resp_all{2},[],1)./sqrt(sum(~isnan(trans_base_resp_all{2}),1)), 'og');
xlim([-10 100])
hline(0, '--k')
title(['All tuned- n = ' num2str(sum(~isnan(nanmean(trans_target_resp_all,2)),1))])
suptitle({titleStr, 'Target response tuning by preference'})
print([fnout 'target_resp_tuning' datasetStr '.pdf'], '-dpdf')

figure;
for iOri = 1:4
    subplot(3,2,iOri)
    for i = 1:2
        tempR = trans_target_resp_SM{i,iOri};
        trans_target_resp_SMcrop{i,iOri} = tempR(:,2:end);
    end
    errorbar(base_dirs(2:end), nanmean(trans_target_resp_SMcrop{1,iOri},1), nanstd(trans_target_resp_SMcrop{1,iOri},[],1)./sqrt(sum(~isnan(trans_target_resp_SMcrop{1,iOri}),1)), '-ok');
    hold on
    errorbar(base_dirs(2:end), nanmean(trans_target_resp_SMcrop{2,iOri},1), nanstd(trans_target_resp_SMcrop{2,iOri},[],1)./sqrt(sum(~isnan(trans_target_resp_SMcrop{2,iOri}),1)), '-or');
    title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp_SMcrop{1,iOri},2)),1))])
    hline(0, '--k')
end
for i = 1:2
    tempR = trans_target_resp_all_SM{i};
    trans_target_resp_all_SMcrop{i} = tempR(:,2:end);
    tempR = trans_target_resp_notori_SM{i};
    trans_target_resp_notori_SMcrop{i} = tempR(:,2:end);
end
subplot(3,2,5)
errorbar(base_dirs(2:end), nanmean(trans_target_resp_notori_SMcrop{1},1), nanstd(trans_target_resp_notori_SMcrop{1},[],1)./sqrt(sum(~isnan(trans_target_resp_notori_SMcrop{1}),1)), '-ok');
hold on
errorbar(base_dirs(2:end), nanmean(trans_target_resp_notori_SMcrop{2},1), nanstd(trans_target_resp_notori_SMcrop{2},[],1)./sqrt(sum(~isnan(trans_target_resp_notori_SMcrop{2}),1)), '-or');
hline(0, '--k')
title(['Untuned- n = ' num2str(sum(~isnan(nanmean(trans_target_resp_notori_SMcrop{1},2)),1))])
subplot(3,2,6)
errorbar(base_dirs(2:end), nanmean(trans_target_resp_all_SMcrop{1},1), nanstd(trans_target_resp_all_SMcrop{1},[],1)./sqrt(sum(~isnan(trans_target_resp_all_SMcrop{1}),1)), '-ok');
hold on
errorbar(base_dirs(2:end), nanmean(trans_target_resp_all_SMcrop{2},1), nanstd(trans_target_resp_all_SMcrop{2},[],1)./sqrt(sum(~isnan(trans_target_resp_all_SMcrop{2}),1)), '-or');
hline(0, '--k')
title(['All tuned- n = ' num2str(sum(~isnan(nanmean(trans_target_resp_all_SMcrop{1},2)),1))])
suptitle({titleStr, 'Tuning by outcome- Hits: Black; Misses: Red'})
print([fnout 'target_resp_tuning_HM' datasetStr '.pdf'], '-dpdf')

figure;
for iOri = 1:4
    subplot(3,2,iOri)
    if sum(~isnan(nanmean(trans_target_resp_SMcrop{1,iOri},2)),1) > 4
        hs = cdfplot(nanmean(trans_target_resp_SMcrop{1,iOri},2));
        set(hs, 'Color', 'k')
        hold on
        hm = cdfplot(nanmean(trans_target_resp_SMcrop{2,iOri},2));
        set(hm, 'Color', 'r')
        [h, p] = kstest2(nanmean(trans_target_resp_SMcrop{1,iOri},2), nanmean(trans_target_resp_SMcrop{2,iOri},2));
    else
        p = NaN;
    end
    title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp_SMcrop{1,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    xlim([-0.05 0.05])
end
subplot(3,2,5)
if sum(~isnan(nanmean(trans_target_resp_notori_SMcrop{1},2)),1) > 4
    hs = cdfplot(nanmean(trans_target_resp_notori_SMcrop{1},2));
    set(hs, 'Color', 'k')
    hold on
    hm = cdfplot(nanmean(trans_target_resp_notori_SMcrop{2},2));
    set(hm, 'Color', 'r')
    [h, p] = kstest2(nanmean(trans_target_resp_notori_SMcrop{1},2),nanmean(trans_target_resp_notori_SMcrop{2},2));
else
    p = NaN;
end
title(['Untuned; n = ' num2str(sum(~isnan(nanmean(trans_target_resp_notori_SMcrop{1},2)),1)) '; p = ' num2str(chop(p,2))])
xlim([-0.05 0.05])
subplot(3,2,6)
if sum(~isnan(nanmean(trans_target_resp_all_SMcrop{1},2)),1) > 4
    hs = cdfplot(nanmean(trans_target_resp_all_SMcrop{1},2));
    set(hs, 'Color', 'k')
    hold on
    hm = cdfplot(nanmean(trans_target_resp_all_SMcrop{2},2));
    set(hm, 'Color', 'r')
    [h, p] = kstest2(nanmean(trans_target_resp_all_SMcrop{1},2),nanmean(trans_target_resp_all_SMcrop{2},2));
else
    p = NaN;
end
title(['All tuned; n = ' num2str(sum(~isnan(nanmean(trans_target_resp_all_SMcrop{1},2)),1)) '; p = ' num2str(chop(p,2))])
xlim([-0.05 0.05])
suptitle({titleStr, 'Hits: Black; Misses: Red'})
print([fnout 'cum_resp_bypref_HM' datasetStr '.pdf'], '-dpdf')

figure;
for iOri = 1:4
    subplot(3,2,iOri)
    if sum(~isnan(nanmean(trans_base_resp_FC{1,iOri},2)),1)>4
        hf = cdfplot(trans_base_resp_FC{1,iOri});
        set(hf, 'Color', 'c')
        hold on
        hc = cdfplot(trans_base_resp_FC{2,iOri});
        set(hc, 'Color', 'b')
        [h, p] = kstest2(trans_base_resp_FC{1,iOri}, trans_base_resp_FC{2,iOri});
    else
        p = NaN;
    end
    title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_base_resp_FC{1,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    xlim([-0.05 0.05])
end
subplot(3,2,5)
if sum(~isnan(nanmean(trans_base_resp_notori_FC{1},2)),1)>4
    hs = cdfplot(trans_base_resp_notori_FC{1});
    set(hs, 'Color', 'c')
    hold on
    hm = cdfplot(trans_base_resp_notori_FC{2});
    set(hm, 'Color', 'b')
    [h, p] = kstest2(trans_base_resp_notori_FC{1},trans_base_resp_notori_FC{2});
else
    p = NaN;
end
title(['Untuned; n = ' num2str(sum(~isnan(nanmean(trans_base_resp_notori_FC{1},2)),1)) '; p = ' num2str(chop(p,2))])
xlim([-0.05 0.05])
subplot(3,2,6)
if sum(~isnan(nanmean(trans_base_resp_all_FC{1},2)),1)>4
    hs = cdfplot(trans_base_resp_all_FC{1});
    set(hs, 'Color', 'c')
    hold on
    hm = cdfplot(trans_base_resp_all_FC{2});
    set(hm, 'Color', 'b')
    [h, p] = kstest2(trans_base_resp_all_FC{1},trans_base_resp_all_FC{2});
else
    p = NaN;
end
title(['All tuned; n = ' num2str(sum(~isnan(nanmean(trans_base_resp_all_FC{1},2)),1)) '; p = ' num2str(chop(p,2))])
suptitle({titleStr, 'FAs: Cyan; CRs: Blue'})
print([fnout 'cum_resp_bypref_FC' datasetStr '.pdf'], '-dpdf')


figure;
for iOri = 1:4
    subplot(3,2,iOri)
    if sum(~isnan(nanmean(trans_base_resp_FC{1,iOri},2)),1) > 4
        hf = cdfplot(trans_base_resp_FC{1,iOri});
        set(hf, 'Color', 'c')
        hold on
        targ = trans_target_resp{1,iOri};
        ha = cdfplot(targ(:,1));
        set(ha, 'Color', 'g')
        [h, p] = kstest2(trans_base_resp_FC{1,iOri}, targ(:,1));
    else
        p = NaN;
    end
    title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_base_resp_FC{1,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    xlim([-0.05 0.05])
end
subplot(3,2,5)
if sum(~isnan(nanmean(trans_base_resp_notori_FC{1},2)),1) > 4
    hs = cdfplot(trans_base_resp_notori_FC{1});
    set(hs, 'Color', 'c')
    hold on
    ha = cdfplot(trans_target_resp_notori(:,1));
    set(ha, 'Color', 'g')
    [h, p] = kstest2(trans_base_resp_notori_FC{1},trans_target_resp_notori(:,1));
else
    p = NaN;
end
title(['Untuned cells; n = ' num2str(sum(~isnan(nanmean(trans_base_resp_notori_FC{1},2)),1)) '; p = ' num2str(chop(p,2))])
xlim([-0.05 0.05])
subplot(3,2,6)
if sum(~isnan(nanmean(trans_base_resp_all_FC{1},2)),1) > 4
    hs = cdfplot(trans_base_resp_all_FC{1});
    set(hs, 'Color', 'c')
    hold on
    ha = cdfplot(trans_target_resp_all(:,1));
    set(ha, 'Color', 'g')
    [h, p] = kstest2(trans_base_resp_all_FC{1},trans_target_resp_all(:,1));
else
    p = NaN;
end
title(['All cells; n = ' num2str(sum(~isnan(nanmean(trans_base_resp_all_FC{1},2)),1)) '; p = ' num2str(chop(p,2))])
xlim([-0.05 0.05])
suptitle({titleStr, 'FAs: Cyan; Auditory: Green'})
print([fnout 'cum_resp_bypref_FA' datasetStr '.pdf'], '-dpdf')

figure;
for iOri = 1:4
    subplot(3,2,iOri)
    targ = trans_target_resp{1,iOri};
    if sum(~isnan(nanmean(targ,2)),1)> 4
        hs = cdfplot(nanmean(targ(:,2:end),2));
        set(hs, 'Color', 'k')
        hold on
        ha = cdfplot(targ(:,1));
        set(ha, 'Color', 'g')
        [h, p] = kstest2(nanmean(targ(:,2:end),2), targ(:,1));
    else
        p = NaN;
    end
    title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(targ,2)),1)) '; p = ' num2str(chop(p,2))])
    xlim([-0.05 0.05])
end
subplot(3,2,5)
if sum(~isnan(nanmean(trans_target_resp_notori,2)),1) > 4
    hs = cdfplot(nanmean(trans_target_resp_notori(:,2:end),2));
    set(hs, 'Color', 'k')
    hold on
    ha = cdfplot(trans_target_resp_notori(:,1));
    set(ha, 'Color', 'g')
    [h, p] = kstest2(nanmean(trans_target_resp_notori(:,2:end),2),trans_target_resp_notori(:,1));
else
    p = NaN;
end
title(['Untuned cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp_notori,2)),1)) '; p = ' num2str(chop(p,2))])
xlim([-0.05 0.05])
subplot(3,2,6)
if sum(~isnan(nanmean(trans_target_resp_all,2)),1) > 4
    hs = cdfplot(nanmean(trans_target_resp_all(:,2:end),2));
    set(hs, 'Color', 'k')
    hold on
    ha = cdfplot(trans_target_resp_all(:,1));
    set(ha, 'Color', 'g')
    [h, p] = kstest2(nanmean(trans_target_resp_all(:,2:end),2),trans_target_resp_all(:,1));
else
    p = NaN;
end
title(['All cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp_all(:,1),2)),1)) '; p = ' num2str(chop(p,2))])
xlim([-0.05 0.05])
suptitle({titleStr, 'Visual: Black; Auditory: Green'})
print([fnout 'cum_resp_bypref_AV' datasetStr '.pdf'], '-dpdf')
    
    