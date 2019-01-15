clear all
close all
CRP_expt_list_all
plotAvg = 0;
area_list = {'C1','C2','LS'};
for id = 1:3
    fprintf(['Day: ' num2str(id) '\n'])
    lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';
    nexp = size(expt(id).date,1);
    all_rew = [];
    all_omit = [];
    all_unexp = [];
    all_earlylick_rew = [];
    all_latelick_rew = [];
    all_earlylick_omit = [];
    all_latelick_omit = [];
    all_earlytrial_rew = [];
    all_latetrial_rew = [];
    all_earlytrial_omit = [];
    all_latetrial_omit = [];
    all_earlytrial_unexp = [];
    all_latetrial_unexp = [];
    all_earlytrial_rew_df = [];
    all_latetrial_rew_df = [];
    all_earlytrial_omit_df = [];
    all_latetrial_omit_df = [];
    all_earlytrial_unexp_df = [];
    all_latetrial_unexp_df = [];
    all_earlytrial_rew_lick = [];
    all_latetrial_rew_lick = [];
    all_earlytrial_omit_lick = [];
    all_latetrial_omit_lick = [];
    all_earlytrial_unexp_lick = [];
    all_latetrial_unexp_lick = [];
    all_short_omit = [];
    all_long_omit = [];
    all_short_unexp = [];
    all_long_unexp = [];
    all_preomit = [];
    all_postomit = [];
    all_early_rew_time = zeros(1,nexp);
    all_late_rew_time = zeros(1,nexp);
    all_early_omit_time = zeros(1,nexp);
    all_late_omit_time = zeros(1,nexp);
    all_rew_df = [];
    all_omit_df = [];
    all_unexp_df = [];
    all_lick_rew = [];
    all_lick_omit = [];
    all_lick_unexp = [];
    all_postrew_lick_rew = [];
    all_postrew_lick_omit = [];
    all_postrew_lick_unexp = [];
    all_early_postrew_lick_rew = [];
    all_late_postrew_lick_rew = [];
    all_postrew_lick_rew_lick = [];
    all_postrew_lick_omit_lick = [];
    all_postrew_lick_unexp_lick = [];
    all_early_postrew_lick_rew_lick = [];
    all_late_postrew_lick_rew_lick = [];
    all_precue_burst = [];
    all_precue_single = [];
    all_precue_burst_df = [];
    all_precue_single_df = [];
    all_precue_single_expt = [];
    all_precue_burst_expt = [];
    all_lastLick_preRew_omit = [];
    all_firstLick_postRew_omit = [];
    all_lastLick_preRew_rew = [];
    all_firstLick_postRew_rew = [];
    all_lastLick_preRew_unexp = [];
    all_firstLick_postRew_unexp = [];
    all_firstPostRewLickEvents = [];
    all_expt_bin = [];
    all_firstPostRewLickEvents_omit = [];
    all_expt_omit_bin = [];
    all_firstPostRewLickEvents_unexp = [];
    all_expt_unexp_bin = [];
    mouse_str = [];
    all_area_id = [];
    expt_areas = zeros(length(area_list),nexp);
    expt_rew_peaks = cell(nexp,length(area_list));
    expt_omit_peaks = cell(nexp,length(area_list));
    expt_unexp_peaks = cell(nexp,length(area_list));
    
    for iexp = 1:nexp
        mouse = expt(id).mouse(iexp,:);
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];
        mouse_str = [mouse_str '_' mouse];
        
        load(fullfile(lg_out,img_fn, [img_fn '_targetAlign.mat']))
        load(fullfile(lg_out,img_fn, [img_fn '_cueAlignLick.mat']))
        load(fullfile(lg_out,img_fn, [img_fn '_lickResp_preCue.mat']))
        if exist(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']))
            load(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']))
            area_temp = zeros(1,size(targetAligndFoverF,2));
            for i = 1:2
                x = find(strcmp(area_list, expt(id).areas{iexp}(i)));
                ind = find(maskCat==i);
                area_temp(1,ind) = x;
                expt_areas(x,iexp) = 1;
            end
        else
            x = find(strcmp(area_list, expt(id).areas{iexp}));
            area_temp = ones(1,size(targetAligndFoverF,2)).*x;
            expt_areas(x,iexp) = 1;
        end
        all_area_id = [all_area_id area_temp];
        
        nTrials = size(targetAligndFoverF,3);
        rew_avg_df = nanmean(targetAligndFoverF(:,:,ind_rew),3);
        omit_avg_df = nanmean(targetAligndFoverF(:,:,ind_omit),3);
        unexp_avg_df = nanmean(targetAligndFoverF(:,:,ind_unexp),3);
        
        rew_avg = nanmean(targetAlign_events(:,:,ind_rew),3);
        omit_avg = nanmean(targetAlign_events(:,:,ind_omit),3);
        unexp_avg = nanmean(targetAlign_events(:,:,ind_unexp),3);

        earlylick_rew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_early_rew)),3);
        latelick_rew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_late_rew)),3);
        earlylick_omit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_early_omit)),3);
        latelick_omit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_late_omit)),3);
        earlytrial_rew_avg = nanmean(targetAlign_events(:,:,ind_rew(find(ind_rew<floor(nTrials./2)))),3);
        latetrial_rew_avg = nanmean(targetAlign_events(:,:,ind_rew(find(ind_rew>floor(nTrials./2)))),3);
        earlytrial_omit_avg = nanmean(targetAlign_events(:,:,ind_omit(find(ind_omit<floor(nTrials./2)))),3);
        latetrial_omit_avg = nanmean(targetAlign_events(:,:,ind_omit(find(ind_omit>floor(nTrials./2)))),3);
        earlytrial_unexp_avg = nanmean(targetAlign_events(:,:,ind_unexp(find(ind_unexp<floor(nTrials./2)))),3);
        latetrial_unexp_avg = nanmean(targetAlign_events(:,:,ind_unexp(find(ind_unexp>floor(nTrials./2)))),3);
        
        earlytrial_rew_avg_df = nanmean(targetAligndFoverF(:,:,ind_rew(find(ind_rew<floor(nTrials./2)))),3);
        latetrial_rew_avg_df = nanmean(targetAligndFoverF(:,:,ind_rew(find(ind_rew>floor(nTrials./2)))),3);
        earlytrial_omit_avg_df = nanmean(targetAligndFoverF(:,:,ind_omit(find(ind_omit<floor(nTrials./2)))),3);
        latetrial_omit_avg_df = nanmean(targetAligndFoverF(:,:,ind_omit(find(ind_omit>floor(nTrials./2)))),3);
        earlytrial_unexp_avg_df = nanmean(targetAligndFoverF(:,:,ind_unexp(find(ind_unexp<floor(nTrials./2)))),3);
        latetrial_unexp_avg_df = nanmean(targetAligndFoverF(:,:,ind_unexp(find(ind_unexp>floor(nTrials./2)))),3);
        
        earlytrial_rew_avg_lick = nanmean(lickCueAlign(:,ind_rew(find(ind_rew<floor(nTrials./2)))),2);
        latetrial_rew_avg_lick = nanmean(lickCueAlign(:,ind_rew(find(ind_rew>floor(nTrials./2)))),2);
        earlytrial_omit_avg_lick = nanmean(lickCueAlign(:,ind_omit(find(ind_omit<floor(nTrials./2)))),2);
        latetrial_omit_avg_lick = nanmean(lickCueAlign(:,ind_omit(find(ind_omit>floor(nTrials./2)))),2);
        earlytrial_unexp_avg_lick = nanmean(lickCueAlign(:,ind_unexp(find(ind_unexp<floor(nTrials./2)))),2);
        latetrial_unexp_avg_lick = nanmean(lickCueAlign(:,ind_unexp(find(ind_unexp>floor(nTrials./2)))),2);
        
        short_omit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_omit_short)),3);
        long_omit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_omit_long)),3);
        short_unexp_avg = nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_short)),3);
        long_unexp_avg = nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_long)),3);
        preomit_avg = nanmean(targetAlign_events(:,:,ind_rew_preomit),3);
        postomit_avg = nanmean(targetAlign_events(:,:,ind_rew_postomit),3);
        
        lick_rew_avg = nanmean(lickCueAlign(:,ind_rew),2);
        lick_omit_avg = nanmean(lickCueAlign(:,ind_omit),2);
        lick_unexp_avg = nanmean(lickCueAlign(:,ind_unexp),2);

        postrew_lick_rew_avg = nanmean(postRew_lickAlignEvents(:,:,ind_rew),3);
        postrew_lick_rew_lick_avg = nanmean(postRew_lickAlign(:,ind_rew),2);
        if sum(~isnan(squeeze(postRew_lickAlignEvents(1,1,ind_omit))))>4
            postrew_lick_omit_avg = nanmean(postRew_lickAlignEvents(:,:,ind_omit),3);
            postrew_lick_omit_lick_avg = nanmean(postRew_lickAlign(:,ind_omit),2);
        else
            postrew_lick_omit_avg = nan(size(postrew_lick_rew_avg));
            postrew_lick_omit_lick_avg = nan(size(postrew_lick_rew_lick_avg));
        end
        postrew_lick_unexp_avg = nanmean(postRew_lickAlignEvents(:,:,ind_unexp),3);
        postrew_lick_unexp_lick_avg = nanmean(postRew_lickAlign(:,ind_unexp),2);
        
        postrew_early_lick_rew_avg = nanmean(postRew_lickAlignEvents(:,:,ind_rew(ind_prerew_early_rew)),3);
        postrew_late_lick_rew_avg = nanmean(postRew_lickAlignEvents(:,:,ind_rew(ind_prerew_late_rew)),3);
        postrew_early_lick_rew_lick_avg = nanmean(postRew_lickAlign(:,ind_rew(ind_prerew_early_rew)),2);
        postrew_late_lick_rew_lick_avg = nanmean(postRew_lickAlign(:,ind_rew(ind_prerew_late_rew)),2);
        
        n_rew = size(ind_rew,2);
        n_omit = size(ind_omit,2);
        n_unexp = size(ind_unexp,2);
        nIC = size(rew_avg,2);
        
        if size(precue_single_lick,3)>3
            precue_single_avg = mean(precue_single_lick,3);
            precue_single_df_avg = mean(precue_single_lick_df,3);
            precue_single_expt = iexp.*ones(1,nIC);
        else
            precue_single_avg = nan(length(tl_precue),nIC);
            precue_single_df_avg = nan(length(tl_precue),nIC);
            precue_single_expt = nan(1,nIC);
        end
        if size(precue_lick_burst,3)>3
            precue_burst_avg = mean(precue_lick_burst,3);
            precue_burst_df_avg = mean(precue_lick_burst_df,3);
            precue_burst_expt = iexp.*ones(1,nIC);
        else
            precue_burst_avg = nan(length(tl_precue),nIC);
            precue_burst_df_avg = nan(length(tl_precue),nIC);
            precue_burst_expt = nan(1,nIC);
        end
        
        if sum(~isnan(lastPreRew_lickAlignEvents(1,1,ind_rew)))>9
            lastPreRewAvg_rew = nanmean(lastPreRew_lickAlignEvents(:,:,ind_rew),3);
        else
            lastPreRewAvg_rew = nan(size(lastPreRew_lickAlignEvents(:,:,1)));
        end
        if sum(~isnan(firstPostRew_lickAlignEvents(1,1,ind_rew)))>9
            firstPostRewAvg_rew = nanmean(firstPostRew_lickAlignEvents(:,:,ind_rew),3);
        else
            firstPostRewAvg_rew = nan(size(firstPostRew_lickAlignEvents(:,:,1)));
        end
        
        if length(ind_omit)>9
            if sum(~isnan(lastPreRew_lickAlignEvents(1,1,ind_omit)))>9
                lastPreRewAvg_omit = nanmean(lastPreRew_lickAlignEvents(:,:,ind_omit),3);
            else
                lastPreRewAvg_omit = nan(size(lastPreRew_lickAlignEvents(:,:,1)));
            end
            if sum(~isnan(firstPostRew_lickAlignEvents(1,1,ind_omit)))>9
                firstPostRewAvg_omit = nanmean(firstPostRew_lickAlignEvents(:,:,ind_omit),3);
            else
                firstPostRewAvg_omit = nan(size(firstPostRew_lickAlignEvents(:,:,1)));
            end
        else
            lastPreRewAvg_omit = nan(size(lastPreRew_lickAlignEvents(:,:,1)));
            firstPostRewAvg_omit = nan(size(firstPostRew_lickAlignEvents(:,:,1)));
        end
        
        if length(ind_unexp)>9
            if sum(~isnan(lastPreRew_lickAlignEvents(1,1,ind_unexp)))>9
                lastPreRewAvg_unexp = nanmean(lastPreRew_lickAlignEvents(:,:,ind_unexp),3);
            else
                lastPreRewAvg_unexp = nan(size(lastPreRew_lickAlignEvents(:,:,1)));
            end
            if sum(~isnan(firstPostRew_lickAlignEvents(1,1,ind_unexp)))>9
                firstPostRewAvg_unexp = nanmean(firstPostRew_lickAlignEvents(:,:,ind_unexp),3);
            else
                firstPostRewAvg_unexp = nan(size(firstPostRew_lickAlignEvents(:,:,1)));
            end
        else
            lastPreRewAvg_unexp = nan(size(lastPreRew_lickAlignEvents(:,:,1)));
            firstPostRewAvg_unexp = nan(size(firstPostRew_lickAlignEvents(:,:,1)));
        end
        
        
        if ~exist(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id)])
            mkdir(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id)])
        end
        
        rewdelay_frames = round(600./frameRateHz);
        edges = ceil([0:0.25:2].*frameRateHz);
        [n bin] = histc(firstPostRewLickFrame,edges);
        ind_trials = cell(1,length(n));
        firstPostRewLickEvents = nan(size(targetAlign_events,1),size(targetAlign_events,2),length(n));
        firstPostRewLickEvents_omit = nan(size(targetAlign_events,1),size(targetAlign_events,2),length(n));
        firstPostRewLickEvents_unexp = nan(size(targetAlign_events,1),size(targetAlign_events,2),length(n));
        ind_rew_temp = cell(1,length(n));
        ind_omit_temp = cell(1,length(n));
        ind_unexp_temp = cell(1,length(n));
        expt_bin = zeros(1,length(n));
        expt_omit_bin = zeros(1,length(n));
        expt_unexp_bin = zeros(1,length(n));
        for ibin = 1:length(n)
            ind_trials{ibin} = find(bin==ibin);
            ind_rew_temp{ibin} = intersect(ind_rew, ind_trials{ibin});
            if length(ind_rew_temp{ibin})>4
                firstPostRewLickEvents(:,:,ibin) = nanmean(targetAlign_events(:,:,ind_rew_temp{ibin}),3);
                expt_bin(1,ibin) = 1;
            end
            ind_omit_temp{ibin} = intersect(ind_omit, ind_trials{ibin});
            if length(ind_omit_temp{ibin})>2
                firstPostRewLickEvents_omit(:,:,ibin) = nanmean(targetAlign_events(:,:,ind_omit_temp{ibin}),3);
                expt_omit_bin(1,ibin) = 1;
            end
            ind_unexp_temp{ibin} = intersect(ind_unexp, ind_trials{ibin});
            if length(ind_unexp_temp{ibin})>2
                firstPostRewLickEvents_unexp(:,:,ibin) = nanmean(targetAlign_events(:,:,ind_unexp_temp{ibin}),3);
                expt_unexp_bin(1,ibin) = 1;
            end
        end
        
        for i = 1:length(area_list)
            ind = find(area_temp == i);
            if length(ind)>5
                avg_rew = nanmean(rew_avg(:,ind),2);
                avg_omit = nanmean(omit_avg(:,ind),2);
                avg_unexp = nanmean(unexp_avg(:,ind),2);
                base_avg = mean(avg_rew(ceil(prewin_frames./2):prewin_frames,:),1);
                base_avg_omit = mean(avg_omit(ceil(prewin_frames./2):prewin_frames,:),1);
                base_avg_unexp = mean(avg_unexp(ceil(prewin_frames./2):prewin_frames,:),1);
                base_std = std(avg_rew(ceil(prewin_frames./2):prewin_frames,:),[],1);
                base_std_omit = std(avg_omit(ceil(prewin_frames./2):prewin_frames,:),[],1);
                base_std_unexp = std(avg_unexp(ceil(prewin_frames./2):prewin_frames,:),[],1);
                if sum(find(avg_rew(prewin_frames+1:prewin_frames+(1.5.*frameRateHz),:)>base_avg+(base_std.*4)))
                    [rew_peak_vals expt_rew_peaks{iexp,i}] = findpeaks(avg_rew(prewin_frames+1:prewin_frames+(1.5.*frameRateHz),:),'MinPeakDistance',5, 'MinPeakHeight', base_avg+(base_std.*4), 'MinPeakProminence', base_std.*4);
                end
                if sum(find(avg_omit(prewin_frames+1:prewin_frames+(1.5.*frameRateHz),:)>base_avg_omit+(base_std_omit.*4)))
                    [omit_peak_vals expt_omit_peaks{iexp,i}] = findpeaks(avg_omit(prewin_frames+1:prewin_frames+(1.5.*frameRateHz),:),'MinPeakDistance',5, 'MinPeakHeight', base_avg_omit+(base_std_omit.*4), 'MinPeakProminence', base_std.*4);
                end
                if sum(find(avg_unexp(prewin_frames+1:prewin_frames+(1.5.*frameRateHz),:)>base_avg_unexp+(base_std_unexp.*4)))
                	[unexp_peak_vals expt_unexp_peaks{iexp,i}] = findpeaks(avg_unexp(prewin_frames+1:prewin_frames+(1.5.*frameRateHz),:),'MinPeakDistance',5, 'MinPeakHeight', base_avg_unexp+(base_std_unexp.*4), 'MinPeakProminence', base_std.*4);
                end
                
                if plotAvg
                    figure;
                    for ibin = 1:length(n)
                        subplot(3,3,ibin)
                        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,ind,ind_rew_temp{ibin}),3),2).*(1000./frameRateHz),(nanstd(nanmean(targetAlign_events(:,ind,ind_rew_temp{ibin}),3),[],2)./sqrt(nIC)).*(1000./frameRateHz));
                        hold on
                        if ibin<length(n)
                            title([num2str(chop(edges(ibin).*(1000/frameRateHz),3)) '-'  num2str(chop(edges(ibin+1).*(1000/frameRateHz),3)) ' ms- n=' num2str(length(ind_rew_temp{ibin})) ' tr'])
                            vline(([edges(ibin) edges(ibin+1)]+rewdelay_frames).*(1000/frameRateHz))
                        else
                            title(['After ' num2str(chop(edges(ibin).*(1000/frameRateHz),3)) ' ms- n=' num2str(length(ind_rew_temp{ibin})) ' tr'])
                            vline([edges(ibin)+rewdelay_frames].*(1000/frameRateHz))
                        end
                        xlabel('Time from cue')
                        ylabel('Spike rate (Hz)')
                    end
                    suptitle([mouse ' ' date ' Area ' area_list{i} ' : binned first lick after rew'])
                    print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\' img_fn '_Area' area_list{i} '_firstLicksAfterRewByMouse.pdf'],'-dpdf', '-fillpage')
               

                    
%                     
%                     figure;
%                     subplot(3,2,1)
%                     plot(tt, earlytrial_rew_avg_lick.*(1000./frameRateHz),'k');
%                     hold on
%                     plot(tt, latetrial_rew_avg_lick.*(1000./frameRateHz),'b');
%                     xlabel('Time from cue')
%                     ylabel('Lick rate (Hz)')
%                     title('Reward')
%                     if id ~= 3
%                         subplot(3,2,2)
%                         plot(tt, earlytrial_omit_avg_lick.*(1000./frameRateHz),'k');
%                         hold on
%                         plot(tt, latetrial_omit_avg_lick.*(1000./frameRateHz),'b');
%                         xlabel('Time from cue')
%                         ylabel('Lick rate (Hz)')
%                         title('Omission')
%                     end
%                     if id == 3
%                         subplot(3,2,2)
%                         plot(tt, earlytrial_unexp_avg_lick.*(1000./frameRateHz),'k');
%                         hold on
%                         plot(tt, latetrial_unexp_avg_lick.*(1000./frameRateHz),'b');
%                         xlabel('Time from cue')
%                         ylabel('Lick rate (Hz)')
%                         title('Unexpected reward')
%                     end
%                     subplot(3,1,2)
%                     shadedErrorBar(tt, nanmean(earlylick_rew_avg(:,ind),2).*(1000./frameRateHz),(nanstd(earlylick_rew_avg(:,ind),[],2)./sqrt(length(ind))).*(1000./frameRateHz),'k');
%                     hold on
%                     shadedErrorBar(tt, nanmean(latelick_rew_avg(:,ind),2).*(1000./frameRateHz),(nanstd(latelick_rew_avg(:,ind),[],2)./sqrt(length(ind))).*(1000./frameRateHz),'b');                    
%                     xlabel('Time from cue')
%                     ylabel('Firing rate (Hz)')
%                     title('Reward')
%                     subplot(3,1,3)
%                     if id ~= 3
%                         shadedErrorBar(tt, nanmean(earlylick_omit_avg(:,ind),2).*(1000./frameRateHz),(nanstd(earlylick_omit_avg(:,ind),[],2)./sqrt(length(ind))).*(1000./frameRateHz),'k');
%                         hold on
%                         shadedErrorBar(tt, nanmean(latelick_omit_avg(:,ind),2).*(1000./frameRateHz),(nanstd(latelick_omit_avg(:,ind),[],2)./sqrt(length(ind))).*(1000./frameRateHz),'b');                    
%                         xlabel('Time from cue')
%                         ylabel('Firing rate (Hz)')
%                         title('Omission')
%                     end
%                     suptitle([mouse ' ' date ' Area ' area_list{i} ' : early (black) late (blue)'])
%                     print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\' img_fn '_Area' area_list{i} '_lickingByMouse.pdf'],'-dpdf', '-fillpage')
% 
%                     
%                     figure;
%                     subplot(3,2,1)
%                     shadedErrorBar(tt, nanmean(rew_avg(:,ind),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,ind,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'k');
%                     hold on
%                     if ~isempty(expt_rew_peaks{iexp,i})
%                         vline([expt_rew_peaks{iexp,i}-1].*1000/frameRateHz)
%                         title([num2str(chop([expt_rew_peaks{iexp,i}-1].*1000/frameRateHz,3)') ' ms from Cue'])
%                     else
%                         title('No peaks found')
%                     end
%                     ylabel('Firing rate (Hz)')
%                     xlabel('Time from cue (ms)')
%                     ylim([0 6])
%                     subplot(3,2,2)
%                     shadedErrorBar(tt, nanmean(rew_avg_df(:,ind),2).*(1000./frameRateHz), (nanstd(nanmean(targetAligndFoverF(:,ind,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'k');
%                     ylabel('dF/F')
%                     xlabel('Time from cue (ms)')
%                     ylim([-1 6])
%                     if length(ind_omit)>5
%                         subplot(3,2,3)
%                         shadedErrorBar(tt, nanmean(omit_avg(:,ind),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,ind,ind_omit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'r');
%                         if ~isempty(expt_omit_peaks{iexp,i})
%                             vline([expt_omit_peaks{iexp,i}-1].*1000/frameRateHz)
%                             title([num2str(chop([expt_omit_peaks{iexp,i}-1].*1000/frameRateHz,3)') ' ms from Cue'])
%                         else
%                             title('No peaks found')
%                         end
%                         ylabel('Firing rate (Hz)')
%                         xlabel('Time from cue (ms)')
%                         ylim([0 6])
%                         subplot(3,2,4)
%                         shadedErrorBar(tt, nanmean(omit_avg_df(:,ind),2).*(1000./frameRateHz), (nanstd(nanmean(targetAligndFoverF(:,ind,ind_omit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'r');
%                         ylabel('dF/F')
%                         xlabel('Time from cue (ms)')
%                         ylim([-1 6])
%                     end
%                     if length(ind_unexp)>5
%                         subplot(3,2,5)
%                         shadedErrorBar(tt, nanmean(unexp_avg(:,ind),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,ind,ind_unexp),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'g');
%                         if ~isempty(expt_unexp_peaks{iexp,i})
%                             vline([expt_unexp_peaks{iexp,i}-1].*1000/frameRateHz)
%                             title([num2str(chop([expt_unexp_peaks{iexp,i}-1].*1000/frameRateHz,3)') ' ms from Cue'])
%                         else
%                             title('No peaks found')
%                         end
%                         ylabel('Firing rate (Hz)')
%                         xlabel('Time from cue (ms)')
%                         ylim([0 6])
%                         subplot(3,2,6)
%                         shadedErrorBar(tt, nanmean(unexp_avg_df(:,ind),2).*(1000./frameRateHz), (nanstd(nanmean(targetAligndFoverF(:,ind,ind_unexp),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'g');
%                         ylabel('dF/F')
%                         xlabel('Time from cue (ms)')
%                         ylim([-1 6])
%                     end
%                     suptitle([mouse ' ' date ' Area ' area_list{i}])
%                     print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\' img_fn '_Area' area_list{i} '_psth.pdf'],'-dpdf', '-fillpage')
% 

%             figure;  
%             IC_list = [];
%             start = 1;
%             for ic = 1:nIC
%                 if start > 49
%                     suptitle([date ' ' mouse '- Cell#' num2str(IC_list(1)) ':' num2str(IC_list(end)) '; Trials: ' num2str(n_rew) ' Rew, ' num2str(n_omit) ' Omit, '  num2str(n_unexp) ' Unexp']) 
%                     print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\' img_fn '_Cell#' num2str(IC_list(1)) '-' num2str(IC_list(end)) '_psth.pdf'],'-dpdf', '-fillpage')
%                     figure;
%                     start = 1;
%                     IC_list = [];
%                 end
%                 subplot(7,7,start)
%                 hold on
%                 if length(ind_omit>1)
%                     shadedErrorBar(tt, omit_avg(:,ic).*(1000./frameRateHz), (nanstd(targetAlign_events(:,ic,ind_omit),[],3)./sqrt(length(ind_omit))).*(1000./frameRateHz),'r');
%                     hold on
%                 end
%                 if length(ind_unexp>1)
%                     shadedErrorBar(tt, unexp_avg(:,ic).*(1000./frameRateHz), (nanstd(targetAlign_events(:,ic,ind_unexp),[],3)./sqrt(length(ind_unexp))).*(1000./frameRateHz),'g');
%                     hold on
%                     vline(1100, 'b')
%                 end
%                 shadedErrorBar(tt, rew_avg(:,ic).*(1000./frameRateHz), (nanstd(targetAlign_events(:,ic,ind_rew),[],3)./sqrt(length(ind_rew))).*(1000./frameRateHz),'k');
%                 xlim([-500 2000])
%                 vline([0 600], 'b')
%                 start = start+1;
%                 IC_list = [IC_list ic];
%             end
%             suptitle([date ' ' mouse '- Cell#' num2str(IC_list(1)) ':' num2str(IC_list(end)) '; Trials: ' num2str(n_rew) ' Rew, ' num2str(n_omit) ' Omit, '  num2str(n_unexp) ' Unexp']) 
%             print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\' img_fn '_Cell#' num2str(IC_list(1)) '-' num2str(IC_list(end)) '_psth.pdf'],'-dpdf', '-fillpage')
                end
            end
        end

        all_rew = [all_rew rew_avg];
        all_omit = [all_omit omit_avg];
        all_unexp = [all_unexp unexp_avg];
        all_earlylick_rew = [all_earlylick_rew earlylick_rew_avg];
        all_latelick_rew = [all_latelick_rew latelick_rew_avg];
        all_earlylick_omit = [all_earlylick_omit earlylick_omit_avg];
        all_latelick_omit = [all_latelick_omit latelick_omit_avg];
        all_earlytrial_rew = [all_earlytrial_rew earlytrial_rew_avg];
        all_latetrial_rew = [all_latetrial_rew latetrial_rew_avg];
        all_earlytrial_omit = [all_earlytrial_omit earlytrial_omit_avg];
        all_latetrial_omit = [all_latetrial_omit latetrial_omit_avg];
        all_earlytrial_unexp = [all_earlytrial_unexp earlytrial_unexp_avg];
        all_latetrial_unexp = [all_latetrial_unexp latetrial_unexp_avg];
        all_earlytrial_rew_df = [all_earlytrial_rew_df earlytrial_rew_avg_df];
        all_latetrial_rew_df = [all_latetrial_rew_df latetrial_rew_avg_df];
        all_earlytrial_omit_df = [all_earlytrial_omit_df earlytrial_omit_avg_df];
        all_latetrial_omit_df = [all_latetrial_omit_df latetrial_omit_avg_df];
        all_earlytrial_unexp_df = [all_earlytrial_unexp_df earlytrial_unexp_avg_df];
        all_latetrial_unexp_df = [all_latetrial_unexp_df latetrial_unexp_avg_df];
        all_earlytrial_rew_lick = [all_earlytrial_rew_lick earlytrial_rew_avg_lick];
        all_latetrial_rew_lick = [all_latetrial_rew_lick latetrial_rew_avg_lick];
        all_earlytrial_omit_lick = [all_earlytrial_omit_lick earlytrial_omit_avg_lick];
        all_latetrial_omit_lick = [all_latetrial_omit_lick latetrial_omit_avg_lick];
        all_earlytrial_unexp_lick = [all_earlytrial_unexp_lick earlytrial_unexp_avg_lick];
        all_latetrial_unexp_lick = [all_latetrial_unexp_lick latetrial_unexp_avg_lick];
        all_short_omit = [all_short_omit short_omit_avg];
        all_long_omit = [all_long_omit long_omit_avg];
        all_short_unexp = [all_short_unexp short_unexp_avg];
        all_long_unexp = [all_long_unexp long_unexp_avg];
        all_preomit = [all_preomit preomit_avg];
        all_postomit = [all_postomit postomit_avg];
        all_rew_df = [all_rew_df rew_avg_df];
        all_omit_df = [all_omit_df omit_avg_df];
        all_unexp_df = [all_unexp_df unexp_avg_df];
        all_lick_rew = [all_lick_rew lick_rew_avg];
        all_lick_omit = [all_lick_omit lick_omit_avg];
        all_lick_unexp = [all_lick_unexp lick_unexp_avg];
        all_early_rew_time(:,iexp) = early_rew_time;
        all_late_rew_time(:,iexp) = late_rew_time;
        all_early_omit_time(:,iexp) = early_omit_time;
        all_late_omit_time(:,iexp) = late_omit_time;
        all_postrew_lick_rew = [all_postrew_lick_rew postrew_lick_rew_avg];
        all_postrew_lick_omit = [all_postrew_lick_omit postrew_lick_omit_avg];
        all_postrew_lick_unexp = [all_postrew_lick_unexp postrew_lick_unexp_avg];
        all_early_postrew_lick_rew = [all_early_postrew_lick_rew postrew_early_lick_rew_avg];
        all_late_postrew_lick_rew = [all_late_postrew_lick_rew postrew_late_lick_rew_avg];
        all_postrew_lick_rew_lick = [all_postrew_lick_rew_lick postrew_lick_rew_lick_avg];
        all_postrew_lick_omit_lick = [all_postrew_lick_omit_lick postrew_lick_omit_lick_avg];
        all_postrew_lick_unexp_lick = [all_postrew_lick_unexp_lick postrew_lick_unexp_lick_avg];
        all_early_postrew_lick_rew_lick = [all_early_postrew_lick_rew_lick postrew_early_lick_rew_lick_avg];
        all_late_postrew_lick_rew_lick = [all_late_postrew_lick_rew_lick postrew_late_lick_rew_lick_avg];
        all_precue_burst = [all_precue_burst precue_burst_avg];
        all_precue_burst_df = [all_precue_burst_df precue_burst_df_avg];
        all_precue_single = [all_precue_single precue_single_avg];
        all_precue_single_df = [all_precue_single_df precue_single_df_avg];
        all_precue_single_expt = [all_precue_single_expt precue_single_expt];
        all_precue_burst_expt = [all_precue_burst_expt precue_burst_expt];
        all_lastLick_preRew_omit = [all_lastLick_preRew_omit lastPreRewAvg_omit];
        all_firstLick_postRew_omit = [all_firstLick_postRew_omit firstPostRewAvg_omit];
        all_lastLick_preRew_rew = [all_lastLick_preRew_rew lastPreRewAvg_rew];
        all_firstLick_postRew_rew = [all_firstLick_postRew_rew firstPostRewAvg_rew];
        all_lastLick_preRew_unexp = [all_lastLick_preRew_unexp lastPreRewAvg_unexp];
        all_firstLick_postRew_unexp = [all_firstLick_postRew_unexp firstPostRewAvg_unexp];
        all_firstPostRewLickEvents = [all_firstPostRewLickEvents firstPostRewLickEvents];
        all_expt_bin = [all_expt_bin; expt_bin];
        all_firstPostRewLickEvents_omit = [all_firstPostRewLickEvents_omit firstPostRewLickEvents_omit];
        all_expt_omit_bin = [all_expt_omit_bin; expt_omit_bin];
        all_firstPostRewLickEvents_unexp = [all_firstPostRewLickEvents_unexp firstPostRewLickEvents_unexp];
        all_expt_unexp_bin = [all_expt_unexp_bin; expt_unexp_bin];
    end

    totIC = size(all_area_id,2);
    totExp = sum(expt_areas,2);
    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        totIC_area = length(ind);
        figure;
        subplot(3,3,1)
        shadedErrorBar(tt, nanmean(all_rew_df(:,ind),2), (nanstd(all_rew_df(:,ind),[],2)./sqrt(totIC_area)),'k');
        xlabel('Time from cue')
        ylabel('dF/F')
        title('Reward')
        ylim([-.05 0.1])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(3,3,2)
        shadedErrorBar(tt, nanmean(all_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Reward')
        ylim([0 5])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(3,3,3)
        shadedErrorBar(tt', nanmean(all_lick_rew(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), (nanstd(all_lick_rew(:,find(expt_areas(i,:))),[],2)./sqrt(totExp(i))).*(1000./frameRateHz),'k');
        xlabel('Time from cue')
        ylabel('Lick rate (Hz)')
        title('Reward')
        ylim([0 10])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(3,3,4)
        shadedErrorBar(tt, nanmean(all_omit_df(:,ind),2), (nanstd(all_omit_df(:,ind),[],2)./sqrt(totIC_area)),'r');
        xlabel('Time from cue')
        ylabel('dF/F')
        title('Omission')
        ylim([-.05 0.1])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(3,3,5)
        shadedErrorBar(tt, nanmean(all_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'r');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omission')
        ylim([0 5])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(3,3,6)
        shadedErrorBar(tt, nanmean(all_lick_omit(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), (nanstd(all_lick_omit(:,find(expt_areas(i,:))),[],2)./sqrt(totExp(i))).*(1000./frameRateHz),'r');
        xlabel('Time from cue')
        ylabel('Lick rate (Hz)')
        title('Omission')
        ylim([0 10])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(3,3,7)
        shadedErrorBar(tt, nanmean(all_unexp_df(:,ind),2), (nanstd(all_unexp_df(:,ind),[],2)./sqrt(totIC_area)),'g');
        xlabel('Time from cue')
        ylabel('dF/F')
        title('Unexpected Reward')
        ylim([-.05 0.1])
        xlim([-500 2000])
        vline([600], 'b')
        subplot(3,3,8)
        shadedErrorBar(tt, nanmean(all_unexp(:,ind),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'g');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Unexpected Reward')
        ylim([0 5])
        xlim([-500 2000])
        vline([600], 'b')
        subplot(3,3,9)
        shadedErrorBar(tt, nanmean(all_lick_unexp(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), (nanstd(all_lick_unexp(:,find(expt_areas(i,:))),[],2)./sqrt(totExp(i))).*(1000./frameRateHz),'g');
        xlabel('Time from cue')
        ylabel('Lick rate (Hz)')
        title('Unexpected Reward')
        ylim([0 10])
        xlim([-500 2000])
        vline([600], 'b')
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice, ' num2str(totIC_area) ' dendrites']);
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_' mouse_str '.fig'])    
    end
            
    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        totIC_area = length(ind);
        figure;
        subplot(2,3,1)
        shadedErrorBar(tt, nanmean(all_earlytrial_rew_df(:,ind),2), (nanstd(all_earlytrial_rew_df(:,ind),[],2)./sqrt(totIC_area)),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_latetrial_rew_df(:,ind),2), (nanstd(all_latetrial_rew_df(:,ind),[],2)./sqrt(totIC_area)),'b');
        xlabel('Time from cue')
        ylabel('dF/F')
        title('Reward')
        ylim([-.05 0.1])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(2,3,2)
        shadedErrorBar(tt, nanmean(all_earlytrial_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_earlytrial_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_latetrial_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_latetrial_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Reward')
        ylim([0 6])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(2,3,3)
        plot(tt, nanmean(all_earlytrial_rew_lick(:,find(expt_areas(i,:))),2).*(1000./frameRateHz),'k');
        hold on
        plot(tt, nanmean(all_latetrial_rew_lick(:,find(expt_areas(i,:))),2).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Lick rate (Hz)')
        title('Reward')
        ylim([0 10])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(2,3,4)
        if id == 3
            shadedErrorBar(tt, nanmean(all_earlytrial_unexp_df(:,ind),2), (nanstd(all_earlytrial_unexp_df(:,ind),[],2)./sqrt(totIC_area)),'k');
            hold on
            shadedErrorBar(tt, nanmean(all_latetrial_unexp_df(:,ind),2), (nanstd(all_latetrial_unexp_df(:,ind),[],2)./sqrt(totIC_area)),'b');
            title('Unexpected')
        else
            shadedErrorBar(tt, nanmean(all_earlytrial_omit_df(:,ind),2), (nanstd(all_earlytrial_omit_df(:,ind),[],2)./sqrt(totIC_area)),'k');
            hold on
            shadedErrorBar(tt, nanmean(all_latetrial_omit_df(:,ind),2), (nanstd(all_latetrial_omit_df(:,ind),[],2)./sqrt(totIC_area)),'b');
            title('Omission')
        end
        xlabel('Time from cue')
        ylabel('dF/F')
        ylim([-.05 0.1])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(2,3,5)
        if id == 3
            shadedErrorBar(tt, nanmean(all_earlytrial_unexp(:,ind),2).*(1000./frameRateHz), (nanstd(all_earlytrial_unexp(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(all_latetrial_unexp(:,ind),2).*(1000./frameRateHz), (nanstd(all_latetrial_unexp(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
            title('Unexpected')
        else
            shadedErrorBar(tt, nanmean(all_earlytrial_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_earlytrial_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(all_latetrial_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_latetrial_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
            title('Omission')
        end
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([0 6])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(2,3,6)
        if id == 3
            plot(tt, nanmean(all_earlytrial_unexp_lick(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), 'k');
            hold on
            plot(tt, nanmean(all_latetrial_unexp_lick(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), 'b');
            title('Unexpected')
        else
            plot(tt, nanmean(all_earlytrial_omit_lick(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), 'k');
            hold on
            plot(tt, nanmean(all_latetrial_omit_lick(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), 'b');
            title('Omission')
        end
        xlabel('Time from cue')
        ylabel('Lick rate (Hz)')
        ylim([0 10])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice, ' num2str(totIC_area) ' dendrites- early (black) vs late (blue) trials']);
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_earlyVlateTrial_' mouse_str '.fig'])    
    end

    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        totIC_area = length(ind);
        figure;
        subplot(2,2,1)
        shadedErrorBar(tt, nanmean(all_earlylick_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_earlylick_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_latelick_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_latelick_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        hold on
        errorbarxy(nanmean(all_early_rew_time(:,find(expt_areas(i,:))),2), 0, nanstd(all_early_rew_time(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2)),0, {'ok', 'k', 'k'});
        hold on
        errorbarxy(nanmean(all_late_rew_time(:,find(expt_areas(i,:))),2), 0, nanstd(all_late_rew_time(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2)),0, {'ob', 'b', 'b'});
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Rewarded trials')
        ylim([-1 6])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(2,2,2)
        shadedErrorBar(tt, nanmean(all_earlylick_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_earlylick_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_latelick_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_latelick_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        hold on
        errorbarxy(nanmean(all_early_omit_time(:,find(expt_areas(i,:))),2),0, nanstd(all_early_omit_time(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2)),0, {'ok', 'k', 'k'});
        hold on
        errorbarxy(nanmean(all_late_omit_time(:,find(expt_areas(i,:))),2), 0, nanstd(all_late_omit_time(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2)),0, {'ob', 'b', 'b'});
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omission trials')
        ylim([-1 6])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice, ' num2str(totIC_area) ' dendrites- Early (black) vs late (blue) lick bursts']);
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_EarlyVsLateLicks_' mouse_str '.fig'])     
    end
    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        totIC_area = length(ind);
        figure;
        if size(all_short_omit(:,ind),2)
            subplot(3,1,1)
            shadedErrorBar(tt, nanmean(all_short_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_short_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(all_long_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_long_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([-1 6])
            xlim([-500 2000])
            vline([600], 'b')
            if id == 4
                vline([1100], 'c')
            end
            title('Omits after short (black) and long (blue) interval')
        end
        if size(all_short_unexp(:,ind),2)
            subplot(3,1,2)
            shadedErrorBar(tt, nanmean(all_short_unexp(:,ind),2).*(1000./frameRateHz), (nanstd(all_short_unexp(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(all_long_unexp(:,ind),2).*(1000./frameRateHz), (nanstd(all_long_unexp(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([-1 6])
            xlim([-500 2000])
            vline([600], 'b')
            if id == 4
                vline([1100], 'c')
            end
            title('Unexpected reward after short (black) and long (blue) interval')
        end
        if size(all_preomit(:,ind),2)
            subplot(3,1,3)
            shadedErrorBar(tt, nanmean(all_preomit(:,ind),2).*(1000./frameRateHz), (nanstd(all_preomit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(all_postomit(:,ind),2).*(1000./frameRateHz), (nanstd(all_postomit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([-1 6])
            xlim([-500 2000])
            vline([600], 'b')
            if id == 4
                vline([1100], 'c')
            end
            title('Reward before (black) and after (blue) omit trial')
        end
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice, ' num2str(totIC_area) ' dendrites']);
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_trialByTrialAnalysis' mouse_str '.fig'])
    end
%     %Single cell Latency analysis    
    base_avg = mean(all_rew(1:prewin_frames,:),1);
    base_avg_omit = mean(all_omit(1:prewin_frames,:),1);
    base_avg_unexp = mean(all_unexp(1:prewin_frames,:),1);
    base_std = std(all_rew(1:prewin_frames,:),[],1);
    base_std_omit = std(all_omit(1:prewin_frames,:),[],1);
    base_std_unexp = std(all_unexp(1:prewin_frames,:),[],1);
%     rew_latency = nan(1,totIC);
%     omit_latency = nan(1,totIC);
%     unexp_latency = nan(1,totIC);
%     
%     for iC = 1:totIC
%         rew_ind = find(all_rew(prewin_frames+1:prewin_frames+postwin_frames,iC)>base_avg(:,iC)+(base_std(:,iC).*1))';
%         omit_ind = find(all_omit(prewin_frames+1:prewin_frames+postwin_frames,iC)>base_avg_omit(:,iC)+(base_std_omit(:,iC).*1))';
%         unexp_ind = find(all_unexp(prewin_frames+1:prewin_frames+postwin_frames,iC)>base_avg_unexp(:,iC)+(base_std_unexp(:,iC).*1))';
%         if strfind(diff(rew_ind), [1 1])
%             temp = strfind(diff(rew_ind), [1 1]);
%             rew_latency(:,iC) = rew_ind(temp(1)).*(1000/frameRateHz);
%         end
%         if strfind(diff(omit_ind), [1 1])
%             temp = strfind(diff(omit_ind), [1 1]);
%             omit_latency(:,iC) = omit_ind(temp(1)).*(1000/frameRateHz);
%         end
%         if strfind(diff(unexp_ind), [1 1])
%             temp = strfind(diff(unexp_ind), [1 1]);
%             unexp_latency(:,iC) = unexp_ind(temp(1)).*(1000/frameRateHz);
%         end
%     end
%     
%     for i = find(sum(expt_areas,2)')
%         ind = find(all_area_id==i);    
%         if id == 1
%             rew_latency_d1 = rew_latency;
%             omit_latency_d1 = omit_latency;
%             unexp_latency_d1 = unexp_latency;
%             ind_d1{i} = find(all_area_id==i);
%             figure;
%             subplot(1,2,1)
%             scatter(rew_latency(:,ind), omit_latency(:,ind), 'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
%             hold on
%             ind_int = intersect(intersect(find(~isnan(rew_latency)), find(~isnan(omit_latency))),ind);
%             n = length(ind_int);
%             errorbarxy(nanmean(rew_latency(:,ind_int),2), nanmean(omit_latency(:,ind_int),2), nanstd(rew_latency(:,ind_int),[],2)./sqrt(n), nanstd(omit_latency(:,ind_int),[],2)./sqrt(n), {'or','r','r'});
%             axis square
%             refline(0,1)
%             xlabel('Reward latency (ms)')
%             ylabel('Omission latency (ms)')
%             title([num2str(n) ' dendrites'])
%             title(['Day ' num2str(id) ' Area ' area_list{i}])
%             savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_RewVOmitLatency' mouse_str '.fig'])
%         end
%         if id == 2
%             figure;
%             subplot(1,2,1)
%             scatter(rew_latency(:,ind), omit_latency(:,ind), 'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
%             hold on
%             ind_int = intersect(intersect(find(~isnan(rew_latency)), find(~isnan(omit_latency))),ind);
%             n = length(ind_int);
%             errorbarxy(nanmean(rew_latency(:,ind_int),2), nanmean(omit_latency(:,ind_int),2), nanstd(rew_latency(:,ind_int),[],2)./sqrt(n), nanstd(omit_latency(:,ind_int),[],2)./sqrt(n), {'or','r','r'});
%             axis square
%             refline(0,1)
%             xlabel('Reward latency (ms)')
%             ylabel('Omission latency (ms)')
%             title([num2str(n) ' dendrites'])
%             suptitle(['Day ' num2str(id) ' Area ' area_list{i}])
%             
%             subplot(1,2,2)
%             errorbar(1:4, [nanmean(rew_latency_d1(:,ind_d1{i}),2) nanmean(omit_latency_d1(:,ind_d1{i}),2) nanmean(rew_latency(:,ind),2)  nanmean(omit_latency(:,ind),2)], [nanstd(rew_latency_d1(:,ind_d1{i}),[],2)./sqrt(sum(~isnan(rew_latency_d1(:,ind_d1{i})),2)) nanstd(omit_latency_d1(:,ind_d1{i}),[],2)./sqrt(sum(~isnan(omit_latency_d1(:,ind_d1{i})),2)) nanstd(rew_latency(:,ind),[],2)./sqrt(sum(~isnan(rew_latency(:,ind)),2)), nanstd(omit_latency(:,ind),[],2)./sqrt(sum(~isnan(omit_latency(:,ind)),2))], 'ok')
%             set(gca, 'Xtick', 1:4, 'Xticklabels',{'RewD1', 'OmitD1', 'RewD2', 'OmitD2'})
%             ylabel('Latency (ms)')
%             xlim([0 5])
%             ylim([0 1500])
%             axis square
%             savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_D1vD2Latency' mouse_str '.fig'])
%         end
%         if id == 3
%             figure;
%             scatter(rew_latency(:,ind), unexp_latency(:,ind), 'ok')
%             hold on
%             ind_int = intersect(intersect(find(~isnan(rew_latency)), find(~isnan(unexp_latency))),ind);
%             n = length(ind_int);
%             errorbarxy(nanmean(rew_latency(:,ind_int),2), nanmean(unexp_latency(:,ind_int),2), nanstd(rew_latency(:,ind_int),[],2)./sqrt(n), nanstd(unexp_latency(:,ind_int),[],2)./sqrt(n), {'og','g','g'});
%             axis square
%             refline(0,1)
%             xlabel('Reward latency (ms)')
%             ylabel('Unexpected latency (ms)')
%             title(['Day ' num2str(id) ' Area ' area_list{i} ': ' num2str(n) ' dendrites'])
%             savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_RewVUnexpLatency' mouse_str '.fig'])
%         end
%     end

        
        
    prerew_resp_cell = [];
    postrew_resp_cell = [];
    postomit_resp_cell = [];
    prerew_supp_cell = [];
    postrew_supp_cell = [];
    postomit_supp_cell = [];
    
    for iC = 1:totIC
        prediff_ind = diff(intersect(prewin_frames+1:prewin_frames+rewdelay_frames, find(all_rew(:,iC)>base_avg(:,iC)+(base_std(:,iC).*1))))';
        if length(prediff_ind)
            if strfind(prediff_ind, [1 1])
                prerew_resp_cell = [prerew_resp_cell iC];
            end
        end
        prediff_ind = diff(intersect(prewin_frames+1:prewin_frames+rewdelay_frames, find(all_rew(:,iC)<base_avg(:,iC)-(base_std(:,iC).*1))))';
        if length(prediff_ind)
            if strfind(prediff_ind, [1 1])
                prerew_supp_cell = [prerew_supp_cell iC];
            end
        end
        postdiff_ind = diff(intersect(prewin_frames+rewdelay_frames+1:prewin_frames+rewdelay_frames+rewdelay_frames, find(all_rew(:,iC)>base_avg(:,iC)+(base_std(:,iC).*1))))';
        if length(postdiff_ind)
            if strfind(postdiff_ind, [1 1])
                postrew_resp_cell = [postrew_resp_cell iC];
            end
        end
        postdiff_ind = diff(intersect(prewin_frames+rewdelay_frames+1:prewin_frames+rewdelay_frames+rewdelay_frames, find(all_rew(:,iC)<base_avg(:,iC)-(base_std(:,iC).*1))))';
        if length(postdiff_ind)
            if strfind(postdiff_ind, [1 1])
                postrew_supp_cell = [postrew_supp_cell iC];
            end
        end
        postomitdiff_ind = diff(intersect(prewin_frames+rewdelay_frames+1:prewin_frames+rewdelay_frames+rewdelay_frames, find(all_omit(:,iC)>base_avg_omit(:,iC)+(base_std_omit(:,iC).*1))))';
        if length(postomitdiff_ind)
            if strfind(postomitdiff_ind, [1 1])
                postomit_resp_cell = [postomit_resp_cell iC];
            end
        end
        postomitdiff_ind = diff(intersect(prewin_frames+rewdelay_frames+1:prewin_frames+rewdelay_frames+rewdelay_frames, find(all_omit(:,iC)<base_avg_omit(:,iC)-(base_std_omit(:,iC).*1))))';
        if length(postomitdiff_ind)
            if strfind(postomitdiff_ind, [1 1])
                postomit_supp_cell = [postomit_supp_cell iC];
            end
        end
    end
    prerew_supp_cell = setdiff(prerew_supp_cell,prerew_resp_cell);
    postrew_supp_cell = setdiff(postrew_supp_cell,postrew_resp_cell);
    postomit_supp_cell = setdiff(postomit_supp_cell,postomit_resp_cell);
    prerew_notresp_cell = setdiff(1:totIC,prerew_resp_cell);
    postrew_notresp_cell = setdiff(1:totIC,postrew_resp_cell);
    postomit_notresp_cell = setdiff(1:totIC,postomit_resp_cell);
    prerew_notsupp_cell = setdiff(1:totIC,prerew_supp_cell);
    postrew_notsupp_cell = setdiff(1:totIC,postrew_supp_cell);
    postomit_notsupp_cell = setdiff(1:totIC,postomit_supp_cell);
   
    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        ind_preresp = intersect(ind, prerew_resp_cell);
        ind_postresp = intersect(ind, postrew_resp_cell);
        ind_notpreresp = intersect(ind, prerew_notresp_cell);
        ind_notpostresp = intersect(ind, postrew_notresp_cell);
        ind_postomitresp = intersect(ind, postomit_resp_cell);
        ind_notpostomitresp = intersect(ind, postomit_notresp_cell);

        totIC_area = length(ind);
        figure; 
        subplot(3,2,1)
        shadedErrorBar(tt, nanmean(all_rew(:,ind_preresp),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind_preresp),[],2)./sqrt(length(ind_preresp))).*(1000./frameRateHz),'k');
        hold on
        if id == 3
            shadedErrorBar(tt, nanmean(all_unexp(:,ind_preresp),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind_preresp),[],2)./sqrt(length(ind_preresp))).*(1000./frameRateHz),'g');
        else
            shadedErrorBar(tt, nanmean(all_omit(:,ind_preresp),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind_preresp),[],2)./sqrt(length(ind_preresp))).*(1000./frameRateHz),'r');
        end
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['Pre-reward responsive: n = ' num2str(length(ind_preresp))])
        ylim([-1 6])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        xlim([-500 2000])
        subplot(3,2,2)
        shadedErrorBar(tt, nanmean(all_rew(:,ind_notpreresp),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind_notpreresp),[],2)./sqrt(length(ind_notpreresp))).*(1000./frameRateHz),'k');
        hold on
        if id == 3
            shadedErrorBar(tt, nanmean(all_unexp(:,ind_notpreresp),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind_notpreresp),[],2)./sqrt(length(ind_notpreresp))).*(1000./frameRateHz),'g');
        else
            shadedErrorBar(tt, nanmean(all_omit(:,ind_notpreresp),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind_notpreresp),[],2)./sqrt(length(ind_notpreresp))).*(1000./frameRateHz),'r');
        end
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['Not pre-reward responsive: n = ' num2str(length(ind_notpreresp))])
        ylim([-1 6])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(3,2,3)
        shadedErrorBar(tt, nanmean(all_rew(:,ind_postresp),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind_postresp),[],2)./sqrt(length(ind_postresp))).*(1000./frameRateHz),'k');
        hold on
        if id == 3
            shadedErrorBar(tt, nanmean(all_unexp(:,ind_postresp),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind_postresp),[],2)./sqrt(length(ind_postresp))).*(1000./frameRateHz),'g');
        else
            shadedErrorBar(tt, nanmean(all_omit(:,ind_postresp),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind_postresp),[],2)./sqrt(length(ind_postresp))).*(1000./frameRateHz),'r');
        end
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['Post-reward responsive: n = ' num2str(length(ind_postresp))])
        ylim([-1 6])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        xlim([-500 2000])
        subplot(3,2,4)
        shadedErrorBar(tt, nanmean(all_rew(:,ind_notpostresp),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind_notpostresp),[],2)./sqrt(length(ind_notpostresp))).*(1000./frameRateHz),'k');
        hold on
        if id == 3
            shadedErrorBar(tt, nanmean(all_unexp(:,ind_notpostresp),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind_notpostresp),[],2)./sqrt(length(ind_notpostresp))).*(1000./frameRateHz),'g');
        else
            shadedErrorBar(tt, nanmean(all_omit(:,ind_notpostresp),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind_notpostresp),[],2)./sqrt(length(ind_notpostresp))).*(1000./frameRateHz),'r');
        end
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['Not post-reward responsive: n = ' num2str(length(ind_notpostresp))])
        ylim([-1 6])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        if id ~=3
            subplot(3,2,5)
            shadedErrorBar(tt, nanmean(all_rew(:,ind_postomitresp),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind_postomitresp),[],2)./sqrt(length(ind_postomitresp))).*(1000./frameRateHz),'k');
            hold on
            if id == 3
                shadedErrorBar(tt, nanmean(all_unexp(:,ind_postomitresp),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind_postomitresp),[],2)./sqrt(length(ind_postomitresp))).*(1000./frameRateHz),'g');
            else
                shadedErrorBar(tt, nanmean(all_omit(:,ind_postomitresp),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind_postomitresp),[],2)./sqrt(length(ind_postomitresp))).*(1000./frameRateHz),'r');
            end
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Post-omit responsive: n = ' num2str(length(ind_postomitresp))])
            ylim([-1 6])
            vline([600], 'b')
            if id == 4
                vline([1100], 'c')
            end
            xlim([-500 2000])
            subplot(3,2,6)
            shadedErrorBar(tt, nanmean(all_rew(:,ind_notpostomitresp),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind_notpostomitresp),[],2)./sqrt(length(ind_notpostomitresp))).*(1000./frameRateHz),'k');
            hold on
            if id == 3
                shadedErrorBar(tt, nanmean(all_unexp(:,ind_notpostomitresp),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind_notpostomitresp),[],2)./sqrt(length(ind_notpostomitresp))).*(1000./frameRateHz),'g');
            else
                shadedErrorBar(tt, nanmean(all_omit(:,ind_notpostomitresp),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind_notpostomitresp),[],2)./sqrt(length(ind_notpostomitresp))).*(1000./frameRateHz),'r');
            end
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Not post-omit responsive: n = ' num2str(length(ind_notpostomitresp))])
            ylim([-1 6])
            xlim([-500 2000])
            vline([600], 'b')
            if id == 4
                vline([1100], 'c')
            end
        end
        if id == 3
            suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice- Reward (black); unexpected (green)'])
        else
            suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice- Reward (black); omit (red)'])
        end
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_preVpostRespCells' mouse_str '.fig'])
    end
   
    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        ind_presupp = intersect(ind, prerew_supp_cell);
        ind_postsupp = intersect(ind, postrew_supp_cell);
        ind_notpresupp = intersect(ind, prerew_notsupp_cell);
        ind_notpostsupp = intersect(ind, postrew_notsupp_cell);
        ind_postomitsupp = intersect(ind, postomit_supp_cell);
        ind_notpostomitsupp = intersect(ind, postomit_notsupp_cell);

        totIC_area = length(ind);
        figure; 
        subplot(3,2,1)
        shadedErrorBar(tt, nanmean(all_rew(:,ind_presupp),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind_presupp),[],2)./sqrt(length(ind_presupp))).*(1000./frameRateHz),'k');
        hold on
        if id == 3
            shadedErrorBar(tt, nanmean(all_unexp(:,ind_presupp),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind_presupp),[],2)./sqrt(length(ind_presupp))).*(1000./frameRateHz),'g');
        else
            shadedErrorBar(tt, nanmean(all_omit(:,ind_presupp),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind_presupp),[],2)./sqrt(length(ind_presupp))).*(1000./frameRateHz),'r');
        end
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['Pre-reward suppressed: n = ' num2str(length(ind_presupp))])
        ylim([-1 6])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        xlim([-500 2000])
        subplot(3,2,2)
        shadedErrorBar(tt, nanmean(all_rew(:,ind_notpresupp),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind_notpresupp),[],2)./sqrt(length(ind_notpresupp))).*(1000./frameRateHz),'k');
        hold on
        if id == 3
            shadedErrorBar(tt, nanmean(all_unexp(:,ind_notpresupp),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind_notpresupp),[],2)./sqrt(length(ind_notpresupp))).*(1000./frameRateHz),'g');
        else
            shadedErrorBar(tt, nanmean(all_omit(:,ind_notpresupp),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind_notpresupp),[],2)./sqrt(length(ind_notpresupp))).*(1000./frameRateHz),'r');
        end
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['Not pre-reward suppressed: n = ' num2str(length(ind_notpresupp))])
        ylim([-1 6])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(3,2,3)
        shadedErrorBar(tt, nanmean(all_rew(:,ind_postsupp),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind_postsupp),[],2)./sqrt(length(ind_postsupp))).*(1000./frameRateHz),'k');
        hold on
        if id == 3
            shadedErrorBar(tt, nanmean(all_unexp(:,ind_postsupp),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind_postsupp),[],2)./sqrt(length(ind_postsupp))).*(1000./frameRateHz),'g');
        else
            shadedErrorBar(tt, nanmean(all_omit(:,ind_postsupp),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind_postsupp),[],2)./sqrt(length(ind_postsupp))).*(1000./frameRateHz),'r');
        end
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['Post-reward suppressed: n = ' num2str(length(ind_postsupp))])
        ylim([-1 6])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        xlim([-500 2000])
        subplot(3,2,4)
        shadedErrorBar(tt, nanmean(all_rew(:,ind_notpostsupp),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind_notpostsupp),[],2)./sqrt(length(ind_notpostsupp))).*(1000./frameRateHz),'k');
        hold on
        if id == 3
            shadedErrorBar(tt, nanmean(all_unexp(:,ind_notpostsupp),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind_notpostsupp),[],2)./sqrt(length(ind_notpostsupp))).*(1000./frameRateHz),'g');
        else
            shadedErrorBar(tt, nanmean(all_omit(:,ind_notpostsupp),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind_notpostsupp),[],2)./sqrt(length(ind_notpostsupp))).*(1000./frameRateHz),'r');
        end
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['Not post-reward suppressed: n = ' num2str(length(ind_notpostsupp))])
        ylim([-1 6])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        if id ~=3
            subplot(3,2,5)
            shadedErrorBar(tt, nanmean(all_rew(:,ind_postomitsupp),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind_postomitsupp),[],2)./sqrt(length(ind_postomitsupp))).*(1000./frameRateHz),'k');
            hold on
            if id == 3
                shadedErrorBar(tt, nanmean(all_unexp(:,ind_postomitsupp),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind_postomitsupp),[],2)./sqrt(length(ind_postomitsupp))).*(1000./frameRateHz),'g');
            else
                shadedErrorBar(tt, nanmean(all_omit(:,ind_postomitsupp),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind_postomitsupp),[],2)./sqrt(length(ind_postomitsupp))).*(1000./frameRateHz),'r');
            end
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Post-omit suppressed: n = ' num2str(length(ind_postomitsupp))])
            ylim([-1 6])
            vline([600], 'b')
            if id == 4
                vline([1100], 'c')
            end
            xlim([-500 2000])
            subplot(3,2,6)
            shadedErrorBar(tt, nanmean(all_rew(:,ind_notpostomitsupp),2).*(1000./frameRateHz), (nanstd(all_rew(:,ind_notpostomitsupp),[],2)./sqrt(length(ind_notpostomitsupp))).*(1000./frameRateHz),'k');
            hold on
            if id == 3
                shadedErrorBar(tt, nanmean(all_unexp(:,ind_notpostomitsupp),2).*(1000./frameRateHz), (nanstd(all_unexp(:,ind_notpostomitsupp),[],2)./sqrt(length(ind_notpostomitsupp))).*(1000./frameRateHz),'g');
            else
                shadedErrorBar(tt, nanmean(all_omit(:,ind_notpostomitsupp),2).*(1000./frameRateHz), (nanstd(all_omit(:,ind_notpostomitsupp),[],2)./sqrt(length(ind_notpostomitsupp))).*(1000./frameRateHz),'r');
            end
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Not post-omit suppressed: n = ' num2str(length(ind_notpostomitsupp))])
            ylim([-1 6])
            xlim([-500 2000])
            vline([600], 'b')
            if id == 4
                vline([1100], 'c')
            end
        end
        if id == 3
            suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice- Reward (black); unexpected (green)'])
        else
            suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice- Reward (black); omit (red)'])
        end
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_preVpostSuppCells' mouse_str '.fig'])
    end

    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        totIC_area = length(ind);
        avg_rew = mean(all_rew(:,ind),2);
        [prerew_max prerew_max_ind] = max(avg_rew(prewin_frames+1:prewin_frames+rewdelay_frames+1,:),[],1);
        [postrew_max postrew_max_ind] = max(avg_rew(prewin_frames+rewdelay_frames+1:prewin_frames+rewdelay_frames+1+rewdelay_frames,:),[],1);
        if id~=3
            figure;
            subplot(2,1,1)
            shadedErrorBar(tt, mean(all_rew(:,ind),2).*1000/frameRateHz, (std(all_rew(:,ind),[],2)./sqrt(length(ind))).*1000/frameRateHz,'k');
            hold on
            shadedErrorBar(tt, mean(all_omit(:,ind),2).*1000/frameRateHz, (std(all_omit(:,ind),[],2)./sqrt(length(ind))).*1000/frameRateHz,'r');            
            vline([prerew_max_ind-1 postrew_max_ind+rewdelay_frames-1].*1000/frameRateHz)
            all_rew_sub = (all_rew-mean(all_rew(1:prewin_frames,:),1)).*(1000/frameRateHz);
            all_omit_sub = (all_omit-mean(all_omit(1:prewin_frames,:),1)).*(1000/frameRateHz);
            subplot(2,2,3)
            scatter(mean(all_rew_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1), mean(all_omit_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1))
            hold on
            errorbarxy(mean(mean(all_rew_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),2),mean(mean(all_omit_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),2),std(mean(all_rew_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),[],2)./sqrt(totIC_area), std(mean(all_omit_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),[],2)./sqrt(totIC_area),{'or', 'r','r'});
            ylabel('Omission- FR (Hz)')
            xlabel('Reward- FR (Hz)')
            xlim([-3 10])
            ylim([-3 10])
            axis square
            refline(1,0)
            title('Pre reward window')
            subplot(2,2,4)
            scatter(mean(all_rew_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1), mean(all_omit_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1))
            hold on
            errorbarxy(mean(mean(all_rew_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),2),mean(mean(all_omit_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),2),std(mean(all_rew_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),[],2)./sqrt(totIC_area), std(mean(all_omit_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),[],2)./sqrt(totIC_area),{'or', 'r','r'});
            ylabel('Omission- FR (Hz)')
            xlabel('Reward- FR (Hz)')
            xlim([-3 10])
            ylim([-3 10])
            refline(1,0)
            axis square
            title('Post reward window')
        end
        if id==3
            figure;
            subplot(2,1,1)
            shadedErrorBar(tt, mean(all_rew(:,ind),2).*1000/frameRateHz, (std(all_rew(:,ind),[],2)./sqrt(length(ind))).*1000/frameRateHz,'k');
            hold on
            shadedErrorBar(tt, mean(all_unexp(:,ind),2).*1000/frameRateHz, (std(all_unexp(:,ind),[],2)./sqrt(length(ind))).*1000/frameRateHz,'g');            
            vline([prerew_max_ind-1 postrew_max_ind+rewdelay_frames-1].*1000/frameRateHz)
            subplot(2,2,3)
            all_rew_sub = (all_rew-mean(all_rew(1:prewin_frames,:),1)).*(1000/frameRateHz);
            all_unexp_sub = (all_unexp-mean(all_unexp(1:prewin_frames,:),1)).*(1000/frameRateHz);
            scatter(mean(all_rew_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1), mean(all_unexp_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1))
            hold on
            errorbarxy(mean(mean(all_rew_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),2),mean(mean(all_unexp_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),2),std(mean(all_rew_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),[],2)./sqrt(totIC_area), std(mean(all_unexp_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),[],2)./sqrt(totIC_area),{'og', 'g','g'});
            ylabel('Unexpected Reward- FR (Hz)')
            xlabel('Reward- FR (Hz)')
            xlim([-3 10])
            ylim([-3 10])
            axis square
            refline(1,0)
            title('Pre reward window')
            subplot(2,2,4)
            scatter(mean(all_rew_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1), mean(all_unexp_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1))
            hold on
            errorbarxy(mean(mean(all_rew_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),2),mean(mean(all_unexp_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),2),std(mean(all_rew_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),[],2)./sqrt(totIC_area), std(mean(all_unexp_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),[],2)./sqrt(totIC_area),{'og', 'g','g'});
            ylabel('Unexpected Reward- FR (Hz)')
            xlabel('Reward- FR (Hz)')
            xlim([-3 10])
            ylim([-3 10])
            refline(1,0)
            axis square
            title('Post reward window')
        end
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totIC_area) ' dendrites'])
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_peakSpikingResp' mouse_str '.fig'])
    end
    
    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        totIC_area = length(ind);
        avg_rew_df = mean(all_rew_df(:,ind),2);
        [prerew_max prerew_max_ind] = max(avg_rew_df(prewin_frames+1:prewin_frames+rewdelay_frames+1,:),[],1);
        [postrew_max postrew_max_ind] = max(avg_rew_df(prewin_frames+rewdelay_frames+1:prewin_frames+rewdelay_frames+1+rewdelay_frames,:),[],1);
        if id~=3
            figure;
            subplot(2,1,1)
            shadedErrorBar(tt, mean(all_rew_df(:,ind),2), (std(all_rew_df(:,ind),[],2)./sqrt(length(ind))),'k');
            hold on
            shadedErrorBar(tt, mean(all_omit_df(:,ind),2), (std(all_omit_df(:,ind),[],2)./sqrt(length(ind))),'r');            
            vline([prerew_max_ind-1 postrew_max_ind+rewdelay_frames-1].*1000/frameRateHz)
            all_rew_df_sub = (all_rew_df-mean(all_rew_df(ceil(prewin_frames./2):prewin_frames,:),1));
            all_omit_df_sub = (all_omit_df-mean(all_omit_df(1:prewin_frames,:),1));
            subplot(2,2,3)
            scatter(mean(all_rew_df_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1), mean(all_omit_df_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1))
            hold on
            errorbarxy(mean(mean(all_rew_df_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),2),mean(mean(all_omit_df_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),2),std(mean(all_rew_df_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),[],2)./sqrt(totIC_area), std(mean(all_omit_df_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),[],2)./sqrt(totIC_area),{'or', 'r','r'});
            ylabel('Omission- dF/F')
            xlabel('Reward- dF/F')
            xlim([-0.2 0.5])
            ylim([-0.2 0.5])
            axis square
            refline(1,0)
            title('Pre reward window')
            subplot(2,2,4)
            scatter(mean(all_rew_df_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1), mean(all_omit_df_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1))
            hold on
            errorbarxy(mean(mean(all_rew_df_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),2),mean(mean(all_omit_df_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),2),std(mean(all_rew_df_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),[],2)./sqrt(totIC_area), std(mean(all_omit_df_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),[],2)./sqrt(totIC_area),{'or', 'r','r'});
            ylabel('Omission- dF/F')
            xlabel('Reward- dF/F')
            xlim([-0.2 0.5])
            ylim([-0.2 0.5])
            refline(1,0)
            axis square
            title('Post reward window')
        end
        if id==3
            figure;
            subplot(2,1,1)
            shadedErrorBar(tt, mean(all_rew_df(:,ind),2), (std(all_rew_df(:,ind),[],2)./sqrt(length(ind))),'k');
            hold on
            shadedErrorBar(tt, mean(all_unexp_df(:,ind),2), (std(all_unexp_df(:,ind),[],2)./sqrt(length(ind))),'g');            
            vline([prerew_max_ind-1 postrew_max_ind+rewdelay_frames-1].*1000/frameRateHz)
            subplot(2,2,3)
            all_rew_df_sub = (all_rew_df-mean(all_rew_df(1:prewin_frames,:),1));
            all_unexp_df_sub = (all_unexp_df-mean(all_unexp_df(1:prewin_frames,:),1));
            scatter(mean(all_rew_df_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1), mean(all_unexp_df_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1))
            hold on
            errorbarxy(mean(mean(all_rew_df_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),2),mean(mean(all_unexp_df_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),2),std(mean(all_rew_df_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),[],2)./sqrt(totIC_area), std(mean(all_unexp_df_sub(prewin_frames+prerew_max_ind-2:prewin_frames+prerew_max_ind,ind),1),[],2)./sqrt(totIC_area),{'og', 'g','g'});
            ylabel('Unexpected Reward- dF/F')
            xlabel('Reward- dF/F')
            xlim([-0.2 0.5])
            ylim([-0.2 0.5])
            axis square
            refline(1,0)
            title('Pre reward window')
            subplot(2,2,4)
            scatter(mean(all_rew_df_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1), mean(all_unexp_df_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1))
            hold on
            errorbarxy(mean(mean(all_rew_df_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),2),mean(mean(all_unexp_df_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),2),std(mean(all_rew_df_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),[],2)./sqrt(totIC_area), std(mean(all_unexp_df_sub(prewin_frames+rewdelay_frames+postrew_max_ind-2:prewin_frames+rewdelay_frames+postrew_max_ind,ind),1),[],2)./sqrt(totIC_area),{'og', 'g','g'});
            ylabel('Unexpected Reward- dF/F')
            xlabel('Reward- dF/F')
            xlim([-0.2 0.5])
            ylim([-0.2 0.5])
            refline(1,0)
            axis square
            title('Post reward window')
        end
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totIC_area) ' dendrites'])
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_peakDFoverFResp' mouse_str '.fig'])
    end
    
    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        totIC_area = length(ind);
        figure;
        subplot(2,3,1)
        shadedErrorBar(tl, nanmean(all_postrew_lick_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_postrew_lick_rew(:,ind),[],2).*(1000./frameRateHz))./sqrt(size(all_postrew_lick_rew(:,ind),2)),'k');
        title(['Reward (black, n = ' num2str(length(all_postrew_lick_rew(:,ind))) ')'])
        hold on
        if sum(~isnan(all_postrew_lick_omit(1,:)),2)>10
            shadedErrorBar(tl, nanmean(all_postrew_lick_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_postrew_lick_omit(:,ind),[],2).*(1000./frameRateHz))./sqrt(size(all_postrew_lick_omit(:,ind),2)),'r');
            title(['Reward (black, n = ' num2str(length(all_postrew_lick_rew(:,ind))) '); omit (red, n = ' num2str(length(all_postrew_lick_omit(:,ind))) ')'])
        end
        xlabel('Time from lick')
        ylabel('Spike rate (Hz)')
        ylim([0 6])
        
        subplot(2,3,2)
        shadedErrorBar(tl, nanmean(all_early_postrew_lick_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_early_postrew_lick_rew(:,ind),[],2).*(1000./frameRateHz))./sqrt(size(all_early_postrew_lick_rew(:,ind),2)),'k');
        hold on
        shadedErrorBar(tl, nanmean(all_late_postrew_lick_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_late_postrew_lick_rew(:,ind),[],2).*(1000./frameRateHz))./sqrt(size(all_late_postrew_lick_rew(:,ind),2)),'b');
        xlabel('Time from lick')
        ylabel('Spike rate (Hz)')
        ylim([0 6])
        title(['Early licks (black); late licks (blue)'])
        subplot(2,3,4)
        shadedErrorBar(tl, nanmean(all_postrew_lick_rew_lick(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), (nanstd(all_postrew_lick_rew_lick(:,find(expt_areas(i,:))),[],2).*(1000./frameRateHz))./sqrt(size(all_postrew_lick_rew_lick(:,find(expt_areas(i,:))),2)),'k');
        hold on
        if sum(~isnan(all_postrew_lick_omit(1,:)),2)>10
            shadedErrorBar(tl, nanmean(all_postrew_lick_omit_lick(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), (nanstd(all_postrew_lick_omit_lick(:,find(expt_areas(i,:))),[],2).*(1000./frameRateHz))./sqrt(size(all_postrew_lick_omit_lick(:,find(expt_areas(i,:))),2)),'r');
        end
        xlabel('Time from lick')
        ylabel('Lick rate (Hz)')
        ylim([0 inf])
        subplot(2,3,5)
        shadedErrorBar(tl, nanmean(all_early_postrew_lick_rew_lick(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), (nanstd(all_early_postrew_lick_rew_lick(:,find(expt_areas(i,:))),[],2).*(1000./frameRateHz))./sqrt(size(all_early_postrew_lick_rew_lick(:,find(expt_areas(i,:))),2)),'k');
        hold on
        shadedErrorBar(tl, nanmean(all_late_postrew_lick_rew_lick(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), (nanstd(all_late_postrew_lick_rew_lick(:,find(expt_areas(i,:))),[],2).*(1000./frameRateHz))./sqrt(size(all_late_postrew_lick_rew_lick(:,find(expt_areas(i,:))),2)),'b');
        xlabel('Time from lick')
        ylabel('Lick rate (Hz)')
        ylim([0 inf])
        if id == 3
            subplot(2,3,3)
            shadedErrorBar(tl, nanmean(all_postrew_lick_unexp(:,ind),2).*(1000./frameRateHz), (nanstd(all_postrew_lick_unexp(:,ind),[],2).*(1000./frameRateHz))./sqrt(size(all_postrew_lick_unexp(:,ind),2)),'g');
            xlabel('Time from lick')
            ylabel('Spike rate (Hz)')
            ylim([0 6])
            title(['Unexpected reward'])
            subplot(2,3,6)
            shadedErrorBar(tl, nanmean(all_postrew_lick_unexp_lick(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), (nanstd(all_postrew_lick_unexp_lick(:,find(expt_areas(i,:))),[],2).*(1000./frameRateHz))./sqrt(size(all_postrew_lick_unexp_lick(:,find(expt_areas(i,:))),2)),'g');
            xlabel('Time from lick')
            ylabel('Lick rate (Hz)')
            ylim([0 inf])
        end
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice- Post-reward lick burst aligned spiking'])
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_postRewLickAlignSpiking' mouse_str '.fig'])
    end
    rew_w = cell(1,3);
    rew_p = cell(1,3);
     for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        totIC_area = length(ind);
        avg_rew = mean(all_rew(:,ind),2);
        avg_omit = mean(all_omit(:,ind),2);
        avg_unexp = mean(all_unexp(:,ind),2);
        base_avg = mean(avg_rew(ceil(prewin_frames./2):prewin_frames,:),1);
        base_avg_omit = mean(avg_omit(ceil(prewin_frames./2):prewin_frames,:),1);
        base_avg_unexp = mean(avg_unexp(ceil(prewin_frames./2):prewin_frames,:),1);
        base_std = std(avg_rew(ceil(prewin_frames./2):prewin_frames,:),[],1);
        base_std_omit = std(avg_omit(ceil(prewin_frames./2):prewin_frames,:),[],1);
        base_std_unexp = std(avg_unexp(ceil(prewin_frames./2):prewin_frames,:),[],1);
        [rew_peak_vals rew_peaks rew_w{i} rew_p{i}] = findpeaks(avg_rew(prewin_frames+1:prewin_frames+(1.5.*frameRateHz),:),'MinPeakDistance',5, 'MinPeakHeight', base_avg+(base_std.*4), 'MinPeakProminence', (base_std.*4));
        [omit_peak_vals omit_peaks] = findpeaks(avg_omit(prewin_frames+1:prewin_frames+(1.5.*frameRateHz),:),'MinPeakDistance',5, 'MinPeakHeight', base_avg_omit+(base_std_omit.*4), 'MinPeakProminence', (base_std.*4));
        [unexp_peak_vals unexp_peaks] = findpeaks(avg_unexp(prewin_frames+1:prewin_frames+(1.5.*frameRateHz),:),'MinPeakDistance',5, 'MinPeakHeight', base_avg_unexp+(base_std_unexp.*4), 'MinPeakProminence', (base_std.*4));
        figure;
        subplot(3,1,1)
        shadedErrorBar(tt, mean(all_rew(:,ind),2).*1000/frameRateHz, (std(all_rew(:,ind),[],2)./sqrt(length(ind))).*1000/frameRateHz,'k');
        hold on
        if ~isempty(rew_peaks)
            vline([rew_peaks-1].*1000/frameRateHz)
            title([num2str(chop([rew_peaks-1].*1000/frameRateHz,3)') ' ms from Cue'])
        else
            title('No peaks found')
        end
        ylabel('Firing rate (Hz)')
        xlabel('Time from cue (ms)')
        if id ~=3
            subplot(3,1,2)
            shadedErrorBar(tt, mean(all_omit(:,ind),2).*1000/frameRateHz, (std(all_omit(:,ind),[],2)./sqrt(length(ind))).*1000/frameRateHz,'r');            
             hold on
            if ~isempty(omit_peaks)
                vline([omit_peaks-1].*1000/frameRateHz)
                title([num2str(chop([omit_peaks-1].*1000/frameRateHz,3)') ' ms from Cue'])
            else
                title('No peaks found')
            end
            ylabel('Firing rate (Hz)')
            xlabel('Time from cue (ms)')
        end
        if id == 3
            subplot(3,1,3)
            shadedErrorBar(tt, mean(all_unexp(:,ind),2).*1000/frameRateHz, (std(all_unexp(:,ind),[],2)./sqrt(length(ind))).*1000/frameRateHz,'g');            
            hold on
            if ~isempty(unexp_peaks)
                vline([unexp_peaks-1].*1000/frameRateHz)
                title([num2str(chop([unexp_peaks-1].*1000/frameRateHz,3)') ' ms from Cue'])
            else
                title('No peaks found')
            end
            ylabel('Firing rate (Hz)')
            xlabel('Time from cue (ms)')
        end
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice- Peak times'])
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_peakSpikeTimes' mouse_str '.fig'])
     end

    for i = find(sum(expt_areas,2)')
        figure;
        for iexp = 1:nexp
            subplot(3,1,1)
            temp = expt_rew_peaks{iexp,i};
            scatter(temp.*(frameRateHz), iexp.*ones(size(temp)), 'xk')
            xlim([0 2000])
            ylim([0 nexp+1])
            ylabel('Expt')
            xlabel('Time from cue (ms)')
            title('Reward')
            hold on
            if id ~= 3
                subplot(3,1,2)
                temp = expt_omit_peaks{iexp,i};
                scatter(temp.*(frameRateHz), iexp.*ones(size(temp)), 'xr')
                xlim([0 2000])
                ylim([0 nexp+1])
                ylabel('Expt')
                xlabel('Time from cue (ms)')
                title('Omit')
                hold on
            end
            if id==3
                subplot(3,1,3)
                temp = expt_unexp_peaks{iexp,i};
                scatter(temp.*(frameRateHz), iexp.*ones(size(temp)), 'xg')
                xlim([0 2000])
                ylim([0 nexp+1])
                ylabel('Expt')
                xlabel('Time from cue (ms)')
                hold on
                title('Unexpected reward')
            end
        end
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ' - Peak times'])
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_peakSpikeTimesByMouse' mouse_str '.fig'])
    end
            
        
    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        temp_single = all_precue_single(:,ind);
        totIC_single = sum(~isnan(temp_single(1,:)),2);
        temp_burst = all_precue_burst(:,ind);
        totIC_burst = sum(~isnan(temp_burst(1,:)),2);
        single_notnan_ind = find(~isnan(all_precue_single_expt(:,ind)));
        burst_notnan_ind = find(~isnan(all_precue_burst_expt(:,ind)));
        nexp_single = length(unique(all_precue_single_expt(:,ind(single_notnan_ind))));
        nexp_burst = length(unique(all_precue_burst_expt(:,ind(burst_notnan_ind))));
        figure;
        subplot(2,2,1)
        shadedErrorBar(tl_precue, nanmean(all_precue_single_df(:,ind),2), (nanstd(all_precue_single_df(:,ind),[],2)./sqrt(totIC_single)),'k');
        xlabel('Time from lick (ms)')
        ylabel('dF/F')
        ylim([-0.02 0.1])
        title(['Single- n= ' num2str(nexp_single) ' mice; ' num2str(totIC_single) 'cells'])
        subplot(2,2,2)
        shadedErrorBar(tl_precue, nanmean(all_precue_single(:,ind),2).*(1000./frameRateHz), (nanstd(all_precue_single(:,ind),[],2)./sqrt(totIC_single)).*(1000./frameRateHz),'k');
        xlabel('Time from lick (ms)')
        ylabel('Spike rate (Hz)')
        title('Single')
        ylim([0 6])
        subplot(2,2,3)
        shadedErrorBar(tl_precue, nanmean(all_precue_burst_df(:,ind),2), (nanstd(all_precue_burst_df(:,ind),[],2)./sqrt(totIC_burst)),'k');
        xlabel('Time from lick (ms)')
        ylabel('dF/F')
        ylim([-0.02 0.1])
        title(['Burst- n= ' num2str(nexp_burst) ' mice; ' num2str(totIC_burst) 'cells'])
        subplot(2,2,4)
        shadedErrorBar(tl_precue, nanmean(all_precue_burst(:,ind),2).*(1000./frameRateHz), (nanstd(all_precue_burst(:,ind),[],2)./sqrt(totIC_burst)).*(1000./frameRateHz),'k');
        xlabel('Time from lick (ms)')
        ylabel('Spike rate (Hz)')
        title('Burst')
        ylim([0 6])
        suptitle(['Day ' num2str(id) ' Area' area_list{i} '- Pre-cue lick-aligned spiking'])
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_preCueLickAlignSpiking' mouse_str '.fig'])
    end
    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        totIC_area = length(ind);
        figure;
        subplot(2,1,1)
        temp_totIC_rew = sum(~isnan(all_lastLick_preRew_rew(1,ind)));
        temp_totIC_omit = sum(~isnan(all_lastLick_preRew_omit(1,ind)));
        temp_totIC_unexp = sum(~isnan(all_lastLick_preRew_unexp(1,ind)));
        shadedErrorBar(tl_rew, nanmean(all_lastLick_preRew_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_lastLick_preRew_rew(:,ind),[],2).*(1000./frameRateHz)./sqrt(temp_totIC_rew)),'k');
        hold on
        if sum(~isnan(all_lastLick_preRew_omit(1,ind)))>5
            shadedErrorBar(tl_rew, nanmean(all_lastLick_preRew_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_lastLick_preRew_omit(:,ind),[],2).*(1000./frameRateHz)./sqrt(temp_totIC_omit)),'r');
            omit_str = ['; Omit: n = ' num2str(temp_totIC_omit)];
        else
            omit_str = [];
        end
        if sum(~isnan(all_lastLick_preRew_unexp(1,ind)))>5
            shadedErrorBar(tl_rew, nanmean(all_lastLick_preRew_unexp(:,ind),2).*(1000./frameRateHz), (nanstd(all_lastLick_preRew_unexp(:,ind),[],2).*(1000./frameRateHz)./sqrt(temp_totIC_unexp)),'g');
            unexp_str = ['; Unexp: n = ' num2str(temp_totIC_unexp)];
        else
            unexp_str = [];
        end
        xlabel('Time from lick (ms)')
        ylabel('Spike rate (Hz)')
        title(['Last lick before reward- Rew: n =' num2str(temp_totIC_rew) omit_str unexp_str])
        ylim([0 6])
        subplot(2,1,2)
        temp_totIC_rew = sum(~isnan(all_firstLick_postRew_rew(1,ind)));
        temp_totIC_omit = sum(~isnan(all_firstLick_postRew_omit(1,ind)));
        temp_totIC_unexp = sum(~isnan(all_firstLick_postRew_unexp(1,ind)));
        shadedErrorBar(tl_rew, nanmean(all_firstLick_postRew_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_firstLick_postRew_rew(:,ind),[],2).*(1000./frameRateHz)./sqrt(temp_totIC_rew)),'k');
        hold on
        if sum(~isnan(all_firstLick_postRew_omit(1,ind)))>5
            shadedErrorBar(tl_rew, nanmean(all_firstLick_postRew_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_firstLick_postRew_omit(:,ind),[],2).*(1000./frameRateHz)./sqrt(temp_totIC_omit)),'r');
            omit_str = ['; Omit: n = ' num2str(temp_totIC_omit)];
        else
            omit_str = [];
        end
        if sum(~isnan(all_firstLick_postRew_unexp(1,ind)))>5
            shadedErrorBar(tl_rew, nanmean(all_firstLick_postRew_unexp(:,ind),2).*(1000./frameRateHz), (nanstd(all_firstLick_postRew_unexp(:,ind),[],2).*(1000./frameRateHz)./sqrt(temp_totIC_unexp)),'g');
            unexp_str = ['; Unexp: n = ' num2str(temp_totIC_unexp)];
        else
            unexp_str = [];
        end
        xlabel('Time from lick (ms)')
        ylabel('Spike rate (Hz)')
        title(['First lick after reward- Rew: n =' num2str(temp_totIC_rew) omit_str unexp_str])
        ylim([0 6])
        suptitle(['Day ' num2str(id) ' Area' area_list{i} '- Lick-aligned spiking'])
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_preAndPostRewLickAlignSpiking' mouse_str '.fig'])
    end
    
    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        figure;
        expt_ind = find(expt_areas(i,:));
        for ibin = 1:length(n)
            subplot(3,3,ibin)
            shadedErrorBar(tt, nanmean(all_firstPostRewLickEvents(:,ind,ibin),2).*(1000./frameRateHz),(nanstd(all_firstPostRewLickEvents(:,ind,ibin),[],2)./sqrt(nIC)).*(1000./frameRateHz));
            hold on
            if ibin<length(n)
                title([num2str(chop(edges(ibin).*(1000/frameRateHz),3)) '-'  num2str(chop(edges(ibin+1).*(1000/frameRateHz),3)) ' ms- n=' num2str(sum(all_expt_bin(expt_ind,ibin),1)) ' mice'])
                vline(([edges(ibin) edges(ibin+1)]+rewdelay_frames).*(1000/frameRateHz))
            else
                title(['After ' num2str(chop(edges(ibin).*(1000/frameRateHz),3)) ' ms- n=' num2str(sum(all_expt_bin(expt_ind,ibin),1)) ' mice'])
                vline(([edges(ibin)]+rewdelay_frames).*(1000/frameRateHz))
            end
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
        end
        suptitle(['Day' num2str(id) ' Area ' area_list{i} ': Reward trials- binned first lick after rew'])
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_binnedFirstLickPostRewSpiking_Rew' mouse_str '.fig'])
        if id ~=3
            figure;
            for ibin = 1:length(n)
                subplot(3,3,ibin)
                shadedErrorBar(tt, nanmean(all_firstPostRewLickEvents_omit(:,ind,ibin),2).*(1000./frameRateHz),(nanstd(all_firstPostRewLickEvents_omit(:,ind,ibin),[],2)./sqrt(nIC)).*(1000./frameRateHz),'r');
                hold on
                if ibin<length(n)
                    title([num2str(chop(edges(ibin).*(1000/frameRateHz),3)) '-'  num2str(chop(edges(ibin+1).*(1000/frameRateHz),3)) ' ms- n=' num2str(sum(all_expt_omit_bin(expt_ind,ibin),1)) ' mice'])
                    vline(([edges(ibin) edges(ibin+1)]+rewdelay_frames).*(1000/frameRateHz))
                else
                    title(['After ' num2str(chop(edges(ibin).*(1000/frameRateHz),3)) ' ms- n=' num2str(sum(all_expt_omit_bin(expt_ind,ibin),1)) ' mice'])
                    vline(([edges(ibin)]+rewdelay_frames).*(1000/frameRateHz))
                end
                xlabel('Time from cue')
                ylabel('Spike rate (Hz)')
            end
            suptitle(['Day' num2str(id) ' Area ' area_list{i} ': Omit trials- binned first lick after rew'])
            savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_binnedFirstLickPostRewSpiking_Omit' mouse_str '.fig'])
        end
        if id == 3
            figure;
            for ibin = 1:length(n)
                subplot(3,3,ibin)
                shadedErrorBar(tt, nanmean(all_firstPostRewLickEvents_unexp(:,ind,ibin),2).*(1000./frameRateHz),(nanstd(all_firstPostRewLickEvents_unexp(:,ind,ibin),[],2)./sqrt(nIC)).*(1000./frameRateHz),'g');
                hold on
                if ibin<length(n)
                    title([num2str(chop(edges(ibin).*(1000/frameRateHz),3)) '-'  num2str(chop(edges(ibin+1).*(1000/frameRateHz),3)) ' ms- n=' num2str(sum(all_expt_unexp_bin(expt_ind,ibin),1)) ' mice'])
                    vline(([edges(ibin) edges(ibin+1)]+rewdelay_frames).*(1000/frameRateHz))
                else
                    title(['After ' num2str(chop(edges(ibin).*(1000/frameRateHz),3)) ' ms- n=' num2str(sum(all_expt_unexp_bin(expt_ind,ibin),1)) ' mice'])
                    vline(([edges(ibin)]+rewdelay_frames).*(1000/frameRateHz))
                end
                xlabel('Time from cue')
                ylabel('Spike rate (Hz)')
            end
            suptitle(['Day' num2str(id) ' Area ' area_list{i} ': Unexpected reward- binned first lick after rew'])
            savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_Summary_binnedFirstLickPostRewSpiking_Unexp' mouse_str '.fig'])
        end
    end


       
        
            
end