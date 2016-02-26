function plotAVTargetSummary(datasetStr)
%takes mouse structure and plots transient responses to target and base
%stimuli
close all
eval(['awFSAVdatasets' datasetStr])
titleStr = datasetStr;
if strcmp(titleStr, '')
    titleStr = 'V1';
else
    titleStr = titleStr(2:end);
end
rc = behavConstsAV;
if strcmp(rc.name,'ashley')
    dataGroup = ['awFSAVdatasets' datasetStr];
else
    dataGroup = [];
end
av = behavParamsAV;
% str = unique({expt.SubNum});
% values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
% mouse_str = [];
% for imouse = 1:size(str,2)
%     mouse_str = [mouse_str 'i' str{imouse} '_'];  
% end
str = unique({expt.SubNum});
values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
mouse_str = ['i' strjoin(str,'_i')];
mouse_ind = find(intersect(cell2mat({av.mouse}),values));
load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' datasetStr '.mat']));
pre_win = mouse(1).expt(1).win(1).frames;
trans_win = mouse(1).expt(1).win(2).frames;
pre_event_frames = mouse(1).expt(1).pre_event_frames;
post_event_frames = mouse(1).expt(1).post_event_frames;
titleStr = [titleStr mouse(1).expt(1).cells(1).name]; 
fnout = fullfile(rc.caOutputDir, dataGroup, [titleStr '_' date '_' mouse_str]); %% maybe lose mouse_str

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


set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
% set(0,'DefaultaxesFontSize', 16)
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
        subplot(n,n2,i)
        i = i+1;
        good_ind = mouse(imouse).expt(iexp).cells(1).ind;
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(3).resp(:,good_ind,:),2),3), 'c');
        hold on
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).resp(:,good_ind,:),2),3), 'g');
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).resp(:,good_ind,:),2),3), 'k');
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
subplot(n,n2,i)
plot(tt, nanmean(resp,2), 'g')
hold on
plot(tt, nanmean(resp_aud,2), 'k')
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
        subplot(n,n2,i)
        i = i+1;
        good_ind = mouse(imouse).expt(iexp).cells(1).ind;
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).resp(:,good_ind,:),2),3), 'g');
        hold on
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).resp(:,good_ind,:),2),3), 'k');
        vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
        vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
        xlim([-10 30])
        resp = cat(2,resp, mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).resp(:,good_ind,:),3));
        resp_aud = cat(2,resp_aud, mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).resp(:,good_ind,:),3));
        title({mouse(imouse).expt(iexp).date, ['n = ' num2str(length(good_ind)) ' cells']})
    end
end
subplot(n,n2,i)
plot(tt,nanmean(resp,2), 'g')
hold on
plot(tt,nanmean(resp_aud,2), 'k')
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
xlim([-10 30])
title(['All cells; n = ' num2str(size(resp,2))])
suptitle({titleStr, 'Visual: Green; Auditory: Black'})
print([fnout 'press_align_TCs' datasetStr '.pdf'], '-dpdf')


for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        for ialign = 2
            for iav = 1:2
                for io = [1 2 5 6]
                    for iDir = 1:size(mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).stimResp,2)
                        tempresp = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).stimResp{iDir};
                        trans_resp = mean(tempresp(trans_win,:,:),1)-mean(tempresp(pre_win,:,:),1);
                        mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).trans_resp(:,iDir) =  squeeze(mean(trans_resp,3));
                    end
                end
                for iDir = 1:size(mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(1).stimResp,2)
                    tempresp = cat(3, mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(1).stimResp{iDir}, mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(2).stimResp{iDir});
                    trans_resp = mean(tempresp(trans_win,:,:),1)-mean(tempresp(pre_win,:,:),1);
                    mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(7).trans_resp(:,iDir) =  squeeze(mean(trans_resp,3));
                    tempresp = cat(3, mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(5).stimResp{iDir}, mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(6).stimResp{iDir});
                    trans_resp = mean(tempresp(trans_win,:,:),1)-mean(tempresp(pre_win,:,:),1);
                    mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(8).trans_resp(:,iDir) =  squeeze(mean(trans_resp,3));
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
                for io = [1 2 5 6]
                	tempresp = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).resp;
                    trans_resp = mean(tempresp(trans_win,:,:),1)-mean(tempresp(pre_win,:,:),1);
                    mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).trans_resp(1,:) =  squeeze(mean(trans_resp,3));
                end
                tempresp = cat(3, mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(1).resp, mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(2).resp);
                trans_resp = mean(tempresp(trans_win,:,:),1)-mean(tempresp(pre_win,:,:),1);
                mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(7).trans_resp(1,:) =  squeeze(mean(trans_resp,3));
                tempresp = cat(3, mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(5).resp, mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(6).resp);
                trans_resp = mean(tempresp(trans_win,:,:),1)-mean(tempresp(pre_win,:,:),1);
                mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(8).trans_resp(1,:) =  squeeze(mean(trans_resp,3));
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


trans_target_resp = cell(8,2,6);
trans_base_resp = cell(8,2,6);
ncells = cell(size(mouse,2),4);
nCyclesOn = cell(6,2);

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
            for i = 1:8
                for iav = 1:2
                    if length(ori_ind)>0
                        if and(iav==1,or(i<3,i>4))
                            tempR = nan(length(ori_ind), length(base_dirs));
                            tempR(:,find(isdir)) = mouse(imouse).expt(iexp).align(2).av(iav).outcome(i).trans_resp(ori_ind,dir_ind);
                            trans_target_resp{i, iav,iOri} = [trans_target_resp{i,iav,iOri}; tempR];
                        elseif and(iav==2,or(i<3,i>4))
                            trans_target_resp{i, iav,iOri} = [trans_target_resp{i,iav,iOri}; nanmean(mouse(imouse).expt(iexp).align(2).av(iav).outcome(i).trans_resp(ori_ind,2:end),2)];
                        elseif or(i>2,i<5)
                            trans_target_resp{i, iav,iOri} = [trans_target_resp{i,iav,iOri}; mouse(imouse).expt(iexp).align(2).av(iav).outcome(i).trans_resp(:,ori_ind)'];
                        end
                        if or(i<3,i==6)
                            trans_base_resp{i, iav,iOri} = [trans_base_resp{i, iav,iOri}; mouse(imouse).expt(iexp).align(1).av(iav).outcome(i).trans_resp(:,ori_ind)'];
                        end
                    elseif and(iexp == 1,imouse == 1)
                        if and(iav==1,or(i<3,i>4))
                            trans_target_resp{i, iav,iOri} = nan(1, length(base_dirs));
                            trans_base_resp{i,iav,iOri} = nan(1, 1);
                        else
                            trans_target_resp{i, iav,iOri} = nan(1, 1);
                            trans_base_resp{i,iav,iOri} = nan(1, 1);
                        end
                    end
                end
            end
        end
        for i = 1:8
            for iav = 1:2
                if length(all_ori_ind)>0
                    iOri = 6;
                    if and(iav==1,or(i<3,i>4))
                        tempR = nan(length(all_ori_ind), length(base_dirs));
                        tempR(:,find(isdir)) = mouse(imouse).expt(iexp).align(2).av(iav).outcome(i).trans_resp(all_ori_ind,dir_ind);
                        trans_target_resp{i, iav,iOri} = [trans_target_resp{i,iav,iOri}; tempR];
                    elseif and(iav==2,or(i<3,i>4))
                        trans_target_resp{i, iav,iOri} = [trans_target_resp{i,iav,iOri}; nanmean(mouse(imouse).expt(iexp).align(2).av(iav).outcome(i).trans_resp(all_ori_ind,2:end),2)];
                    elseif or(i>2,i<5)
                        trans_target_resp{i, iav,iOri} = [trans_target_resp{i,iav,iOri}; mouse(imouse).expt(iexp).align(2).av(iav).outcome(i).trans_resp(:,all_ori_ind)'];
                    end
                    if or(i<3,i==6)
                        trans_base_resp{i, iav,iOri} = [trans_base_resp{i, iav,iOri}; mouse(imouse).expt(iexp).align(1).av(iav).outcome(i).trans_resp(:,all_ori_ind)'];
                    end
                elseif and(iexp == 1,imouse == 1)
                    if and(iav==1,or(i<3,i>4))
                        trans_target_resp{i, iav,iOri} = nan(1, length(base_dirs));
                        trans_base_resp{i,iav,iOri} = nan(1, 1);
                    else
                        trans_target_resp{i, iav,iOri} = nan(1, 1);
                        trans_base_resp{i,iav,iOri} = nan(1, 1);
                    end
                end
            end
        end
        notori_ind = setdiff(good_ind, all_ori_ind);
        for i = 1:8
            for iav = 1:2
                if length(notori_ind)>0
                    iOri = 5;
                    if and(iav==1,or(i<3,i>4))
                        tempR = nan(length(notori_ind), length(base_dirs));
                        tempR(:,find(isdir)) = mouse(imouse).expt(iexp).align(2).av(iav).outcome(i).trans_resp(notori_ind,dir_ind);
                        trans_target_resp{i, iav,iOri} = [trans_target_resp{i,iav,iOri}; tempR];
                    elseif and(iav==2,or(i<3,i>4))
                        trans_target_resp{i, iav,iOri} = [trans_target_resp{i,iav,iOri}; nanmean(mouse(imouse).expt(iexp).align(2).av(iav).outcome(i).trans_resp(notori_ind,2:end),2)];
                    elseif or(i>2,i<5)
                        trans_target_resp{i, iav,iOri} = [trans_target_resp{i,iav,iOri}; mouse(imouse).expt(iexp).align(2).av(iav).outcome(i).trans_resp(:,notori_ind)'];
                    end
                    if or(i<3,i==6)
                        trans_base_resp{i, iav,iOri} = [trans_base_resp{i, iav,iOri}; mouse(imouse).expt(iexp).align(1).av(iav).outcome(i).trans_resp(:,notori_ind)'];
                    end
                elseif and(iexp == 1,imouse == 1)
                    if and(iav==1,or(i<3,i>4))
                        trans_target_resp{i, iav,iOri} = nan(1, length(base_dirs));
                        trans_base_resp{i,iav,iOri} = nan(1, 1);
                    else
                        trans_target_resp{i, iav,iOri} = nan(1, 1);
                        trans_base_resp{i,iav,iOri} = nan(1, 1);
                    end
                end
            end
        end
        for iav = 1:2
            for io = 1:6
                ialign = 2;
                nCyclesOn{io,iav} = [nCyclesOn{io,iav} mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(io).tcyc];
            end
        end
    end
end

cell_set{5} = 'Untuned';
cell_set{6} = 'All Tuned';

figure;
for iOri = 1:6
    subplot(3,2,iOri)
    errorbar(base_dirs, nanmean(trans_target_resp{1,1,iOri},1), nanstd(trans_target_resp{1,1,iOri},[],1)./sqrt(sum(~isnan(trans_target_resp{1,1,iOri}),1)), '-ob');
    hold on
    errorbar(0, nanmean(trans_base_resp{1,1,iOri},1), nanstd(trans_base_resp{1,1,iOri},[],1)./sqrt(sum(~isnan(trans_base_resp{1,1,iOri}),1)), 'og');
    errorbar(0, nanmean(trans_base_resp{1,2,iOri},1), nanstd(trans_base_resp{1,2,iOri},[],1)./sqrt(sum(~isnan(trans_base_resp{1,2,iOri}),1)), 'ok');
    if iOri<5
        title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp{1,1,iOri},2)),1))])
    else
        title([cell_set{iOri} ' cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp{1,1,iOri},2)),1))])
    end
    xlim([-10 100])
    ylim([0 0.03])
end
suptitle([titleStr '- Target response tuning by preference- Success only'])
print([fnout 'target_resp_tuning_success' datasetStr '.pdf'], '-dpdf')

% g = flipud(gray(length(base_dirs)));
% figure;
% tt = 1-pre_event_frames:post_event_frames;
% for iOri = 1:4
%     subplot(3,2,iOri)
%     for iDir = 1:length(base_dirs)
%         TC_temp = TC_target_resp{iOri};
%         plot(tt, squeeze(nanmean(TC_temp(:,:,iDir),2)), 'Color', g(iDir,:));
%         hold on
%     end
%     vline(0, '--k')
%     vline([trans_win(1) trans_win(end)]-pre_event_frames,'--r')
%     xlim([-20 20])
%     title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp{iOri},2)),1))])
% end
% subplot(3,2,5)
% for iDir = 1:length(base_dirs)
%     plot(tt, squeeze(nanmean(TC_target_resp_notori(:,:,iDir),2)), 'Color', g(iDir,:));
%     hold on
% end
% vline(0, '--k')
% vline([trans_win(1) trans_win(end)]-pre_event_frames,'--r')
% xlim([-20 20])
% title(['Untuned- n = ' num2str(sum(~isnan(nanmean(TC_target_resp_notori(1,:,:),3)),2))])
% subplot(3,2,6)
% for iDir = 1:length(base_dirs)
%     plot(tt, squeeze(nanmean(TC_target_resp_all(:,:,iDir),2)), 'Color', g(iDir,:));
%     hold on
% end
% vline(0, '--k')
% vline([trans_win(1) trans_win(end)]-pre_event_frames,'--r')
% xlim([-20 20])
% title(['All tuned- n = ' num2str(sum(~isnan(nanmean(TC_target_resp_all(1,:,:),3)),2))])

trans_target_resp_crop = cell(size(trans_target_resp));
for iOri = 1:6
    for ii = 1:8
        for i = 1:2
            tempR = trans_target_resp{ii,i,iOri};
            trans_target_resp_crop{ii,i,iOri} = tempR(:,2:end);
        end
    end
end

figure;
for iOri = 1:6
    subplot(3,2,iOri)
    errorbar(base_dirs(2:end), nanmean(trans_target_resp_crop{5, 1,iOri},1), nanstd(trans_target_resp_crop{5,1,iOri},[],1)./sqrt(sum(~isnan(trans_target_resp_crop{5,1,iOri}),1)), '-ok');
    hold on
    errorbar(base_dirs(2:end), nanmean(trans_target_resp_crop{6, 1,iOri},1), nanstd(trans_target_resp_crop{6,1,iOri},[],1)./sqrt(sum(~isnan(trans_target_resp_crop{6,1,iOri}),1)), '-or');
    if iOri<5
        title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp_crop{5,1,iOri},2)),1))])
    else
        title([cell_set{iOri} ' cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp_crop{5,1,iOri},2)),1))])
    end
    xlim([-10 100])
    ylim([-.01 0.05])
end
suptitle([titleStr '- Target response tuning by preference- Matched Hits and Misses'])
print([fnout 'target_resp_tuning_HM' datasetStr '.pdf'], '-dpdf')

figure;
for iOri = 1:6
    subplot(3,2,iOri)
    if sum(~isnan(nanmean(trans_target_resp_crop{5,1,iOri},2)),1) > 4
        hs = cdfplot(nanmean(trans_target_resp_crop{5,1,iOri},2));
        set(hs, 'Color', 'k')
        hold on
        hm = cdfplot(nanmean(trans_target_resp_crop{6,1,iOri},2));
        set(hm, 'Color', 'r')
        [h, p] = kstest2(nanmean(trans_target_resp_crop{5,1,iOri},2), nanmean(trans_target_resp_crop{6,1,iOri},2));
    else
        p = NaN;
    end
    if iOri<5
        title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp_crop{5,1,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    else
        title([cell_set{iOri} ' cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp_crop{5,1,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    end
    xlim([-0.05 0.05])
end
suptitle([titleStr '- Hits: Black- ' num2str(chop(mean(nCyclesOn{5,1},2),2)) ' cyc; Misses: Red ' num2str(chop(mean(nCyclesOn{6,1},2),2)) ' cyc'])
print([fnout 'cum_resp_bypref_HM' datasetStr '.pdf'], '-dpdf')

figure;
for iOri = 1:6
    subplot(3,2,iOri)
    if sum(~isnan(nanmean(trans_target_resp{3,1,iOri},2)),1)>4
        hf = cdfplot(trans_target_resp{3,1,iOri});
        set(hf, 'Color', 'c')
        hold on
        hc = cdfplot(trans_target_resp{4,1,iOri});
        set(hc, 'Color', 'b')
        [h, p] = kstest2(trans_target_resp{3,1,iOri}, trans_target_resp{4,1,iOri});
    else
        p = NaN;
    end
    if iOri<5
        title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp{3,1,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    else
        title([cell_set{iOri} ' cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp{3,1,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    end
    xlim([-0.05 0.05])
end
suptitle([titleStr '- FAs: Cyan ' num2str(chop(mean(nCyclesOn{3,1},2),2)) ' cyc; CRs: Blue ' num2str(chop(mean(nCyclesOn{4,1},2),2)) ' cyc'])
print([fnout 'cum_resp_bypref_FC' datasetStr '.pdf'], '-dpdf')


figure;
for iOri = 1:6
    subplot(3,2,iOri)
    if sum(~isnan(nanmean(trans_target_resp{3,1,iOri},2)),1) > 4
        hf = cdfplot(trans_target_resp{3,1,iOri});
        set(hf, 'Color', 'c')
        hold on
        ha = cdfplot(trans_target_resp{1,2,iOri});
        set(ha, 'Color', 'k')
        [h, p] = kstest2(trans_target_resp{3,1,iOri}, trans_target_resp{1,2,iOri});
    else
        p = NaN;
    end
    if iOri<5
        title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp{3,1,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    else
        title([cell_set{iOri} ' cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp{3,1,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    end
    xlim([-0.05 0.05])
end
suptitle([titleStr '- FAs: Cyan ' num2str(chop(mean(nCyclesOn{3,1},2),2)) ' cyc; Auditory: Black ' num2str(chop(mean(nCyclesOn{1,2},2),2)) ' cyc'])
print([fnout 'cum_resp_bypref_FA' datasetStr '.pdf'], '-dpdf')

%auditory false alarms
figure;
for iOri = 1:6
    subplot(3,2,iOri)
    if sum(~isnan(nanmean(trans_target_resp{3,2,iOri},2)),1) > 4
        hf = cdfplot(trans_target_resp{3,2,iOri});
        set(hf, 'Color', 'c')
        hold on
        ha = cdfplot(trans_target_resp{1,2,iOri});
        set(ha, 'Color', 'k')
        [h, p] = kstest2(trans_target_resp{3,2,iOri},trans_target_resp{1,2,iOri});
    else
        p = NaN;
    end
    if iOri<5
        title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp{3,2,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    else
        title([cell_set{iOri} ' cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp{3,2,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    end
    xlim([-0.05 0.05])
end
suptitle([titleStr '- Auditory FAs: Cyan ' num2str(chop(mean(nCyclesOn{3,2},2),2)) ' cyc; Auditory: Black ' num2str(chop(mean(nCyclesOn{1,2},2),2)) ' cyc'])
print([fnout 'cum_resp_bypref_FAAud' datasetStr '.pdf'], '-dpdf')

figure;
for iOri = 1:6
    subplot(3,2,iOri)
    if sum(~isnan(nanmean(trans_target_resp_crop{1,1,iOri},2)),1)> 4
        temp = trans_target_resp_crop{1,1,iOri};
        hs = cdfplot(nanmean(temp,2));
        set(hs, 'Color', 'g')
        hold on
        ha = cdfplot(trans_target_resp{1,2,iOri});
        set(ha, 'Color', 'k')
        [h, p] = kstest2(nanmean(temp,2), trans_target_resp{1,2,iOri});
    else
        p = NaN;
    end
    if iOri<5
        title([mouse(1).expt(1).cells(iOri+1).name ' deg cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp_crop{1,1,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    else
        title([cell_set{iOri} ' cells; n = ' num2str(sum(~isnan(nanmean(trans_target_resp_crop{1,1,iOri},2)),1)) '; p = ' num2str(chop(p,2))])
    end
    xlim([-0.05 0.05])
end
suptitle([titleStr '- Visual: Green ' num2str(chop(mean(nCyclesOn{1,1},2),2)) ' cyc; Auditory: Black ' num2str(chop(mean(nCyclesOn{1,2},2),2)) ' cyc'])
print([fnout 'cum_resp_bypref_AV_Hits' datasetStr '.pdf'], '-dpdf')
    
    