clear all
close all
CRP_expt_list_all
for id = 1:4
lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';
nexp = size(expt(id).date,1);
fprintf(['Day ' num2str(id) '\n'])
    for iexp = 1:nexp
        mouse = expt(id).mouse(iexp,:);
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];

        load(fullfile(lg_out,img_fn, [img_fn '_input.mat']))
        load(fullfile(lg_out,img_fn, [img_fn '_targetAlign.mat']))
        
        if id == 4
            rewDelay_frames =  round(1.1.*frameRateHz);
        else
            rewDelay_frames =  round(0.6.*frameRateHz);
        end
        
        cTargetOn = celleqel2mat_padded(input.cTargetOn);
        nTrials = size(cTargetOn,2);
        nIC = size(targetAlign_events,2);
        lickCueAlign =  nan(prewin_frames+postwin_frames,nTrials);
        lickCounterVals = cell(1,nTrials);
        nFrames = input.counterValues{nTrials}(end);
        lickDelay_frames =  round(0.1.*frameRateHz);
        lickSearch_frames =  round(0.3.*frameRateHz);
        
        postLick_frames =  round(0.5.*frameRateHz);
        lickSearchRange = prewin_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames;
        lickBurstStart = nan(1,nTrials);
        postRew_lickBurstStart = nan(1,nTrials);
        postRew_lickAlignEvents = nan(2.*postLick_frames,nIC,nTrials);
        postRew_lickAlign = nan(2.*postLick_frames,nTrials);
        lastPreRew_lickAlignEvents = nan(2.*postLick_frames,nIC,nTrials);
        firstPostRew_lickAlignEvents = nan(2.*postLick_frames,nIC,nTrials);
        lastPreRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
        firstPostRew_lickAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
        lastPreRew_lickAlign = nan(3.*postLick_frames,nTrials);
        firstPostRew_lickAlign = nan(3.*postLick_frames,nTrials);
        rewAlignEvents = nan(3.*postLick_frames,nIC,nTrials);
        preRewTrials = [];
        postRewTrials = [];
        postRew_lickSearchRange = prewin_frames+rewDelay_frames+lickDelay_frames:size(lickCueAlign,1)-lickSearch_frames-postLick_frames;
        lastPreRewLickFrame = nan(1,nTrials);
        firstPostRewLickFrame = nan(1,nTrials);
        for itrial = 1:nTrials
            if cTargetOn(itrial)+postwin_frames-1 < nFrames
                lickTimes = input.lickometerTimesUs{itrial};
                counterTimes = input.counterTimesUs{itrial};
                counterVals = input.counterValues{itrial};
                lickCounterVals{itrial} = zeros(size(lickTimes));
                lickTC{itrial} = zeros(size(counterVals));
                for icount = 1:length(counterTimes)-1
                    ind = find(lickTimes>counterTimes(icount) & lickTimes<counterTimes(icount+1));
                    if length(ind)>0
                        lickCounterVals{itrial}(1,ind) = icount;
                    end
                end
                for ival = 1:length(counterTimes)
                    ind = find(lickCounterVals{itrial} == ival);
                    if length(ind)>0
                        lickTC{itrial}(1,ival) = length(ind);
                    end
                end
                if input.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
                    lickCueAlign(:,itrial) = lickTC{itrial}(1,cTargetOn(itrial)-prewin_frames-counterVals(1):cTargetOn(itrial)+postwin_frames-1-counterVals(1))';
                end
                ind = intersect(lickSearchRange,find(lickCueAlign(:,itrial)));
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,itrial),1) >= 3
                        lickBurstStart(:,itrial) = ilick;
                        break
                    end
                end
                ind = intersect(postRew_lickSearchRange,find(lickCueAlign(:,itrial)));
                for i = 1:length(ind)
                    ilick = ind(i);
                    if sum(lickCueAlign(ilick:ilick+lickSearch_frames-1,itrial),1) >= 3
                        postRew_lickBurstStart(:,itrial) = ilick;
                        postRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ilick-postLick_frames:ilick+postLick_frames-1,:,itrial);
                        postRew_lickAlign(:,itrial) = lickCueAlign(ilick-postLick_frames:ilick+postLick_frames-1,itrial);
                        break
                    end
                end
                ind_pre = find(lickCueAlign(prewin_frames:prewin_frames+rewDelay_frames,itrial),1,'last');
                if ~isempty(ind_pre)
                    lastPreRewLickFrame(1,itrial) = ind_pre;
                    preRewTrials = [preRewTrials itrial];
                    lastPreRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,:,itrial);
                    lastPreRew_lickAlign(:,itrial) = lickCueAlign(ind_pre-postLick_frames+prewin_frames:ind_pre+postLick_frames+postLick_frames-1+prewin_frames,itrial);
                end
                ind_post = find(lickCueAlign(prewin_frames+rewDelay_frames:prewin_frames+postwin_frames-postLick_frames-postLick_frames,itrial),1,'first');
                if ~isempty(ind_post)
                    firstPostRewLickFrame(1,itrial) = ind_post;
                    postRewTrials = [postRewTrials itrial];
                    firstPostRew_lickAlignEvents(:,:,itrial) = targetAlign_events(ind_post-postLick_frames+prewin_frames+rewDelay_frames:ind_post+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,:,itrial);
                    firstPostRew_lickAlign(:,itrial) = lickCueAlign(ind_post-postLick_frames+prewin_frames+rewDelay_frames:ind_post+postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,itrial);
                end
                rewAlignEvents(:,:,itrial) =targetAlign_events(-postLick_frames+prewin_frames+rewDelay_frames:postLick_frames+postLick_frames-1+prewin_frames+rewDelay_frames,:,itrial);
            end
        end

        figure;
        subplot(3,1,1)
        shadedErrorBar(tt, nanmean(lickCueAlign(:,ind_rew),2).*(1000./frameRateHz), (nanstd(lickCueAlign(:,ind_rew),[],2)./sqrt(length(ind_rew))).*(1000./frameRateHz));
        hold on
        scatter((lickBurstStart(:,ind_rew)-prewin_frames).*(1000./frameRateHz), 10.*ones(size(lickBurstStart(:,ind_rew))), 'x');
        xlabel('Time from cue')
        ylabel('Lick rate (Hz)')
        title('Reward')
        if length(ind_omit>5)
        subplot(3,1,2)
        shadedErrorBar(tt, nanmean(lickCueAlign(:,ind_omit),2).*(1000./frameRateHz), (nanstd(lickCueAlign(:,ind_omit),[],2)./sqrt(length(ind_omit))).*(1000./frameRateHz),'r');
        hold on
        scatter((lickBurstStart(:,ind_omit)-prewin_frames).*(1000./frameRateHz), 10.*ones(size(lickBurstStart(:,ind_omit))), 'x');
        xlabel('Time from cue')
        ylabel('Lick rate (Hz)')
        title('Omit')
        end
        if length(ind_unexp>5)
        subplot(3,1,3)
        shadedErrorBar(tt, nanmean(lickCueAlign(:,ind_unexp),2).*(1000./frameRateHz), nanstd(lickCueAlign(:,ind_unexp),[],2)./sqrt(length(ind_unexp)).*(1000./frameRateHz),'g');
        hold on
        scatter((lickBurstStart(:,ind_unexp)-prewin_frames).*(1000./frameRateHz), 10.*ones(size(lickBurstStart(:,ind_unexp))), 'x');
        xlabel('Time from cue')
        ylabel('Lick rate (Hz)')
        title('Unexpected reward')
        end    
        suptitle([date ' ' mouse])
        savefig(fullfile(lg_out,img_fn, [img_fn '_cueAlign_lickHz.fig']))

        nIC = size(targetAlign_events,2);
        if sum(~isnan(lickBurstStart(:,ind_rew)))>6
            [sortlick sortlick_ind] = sort(lickBurstStart(:,ind_rew),'ascend');
            nburst = sum(~isnan(lickBurstStart(:,ind_rew)));
            nnan = sum(isnan(lickBurstStart(:,ind_rew)));
            ind_early_rew = sortlick_ind(1:floor(nburst/4));
            ind_late_rew = sortlick_ind(nburst-floor(nburst/4)+1:end-nnan);
            early_rew_time = mean((lickBurstStart(:,ind_rew(ind_early_rew))-prewin_frames).*(1000./frameRateHz),2);
            late_rew_time = mean((lickBurstStart(:,ind_rew(ind_late_rew))-prewin_frames).*(1000./frameRateHz),2);
        else
            ind_early_rew = [];
            ind_late_rew = [];
            early_rew_time = [];
            late_rew_time = [];
        end
        if sum(~isnan(lickBurstStart(:,ind_omit)))>6
            [sortlick sortlick_ind] = sort(lickBurstStart(:,ind_omit),'ascend');
            nburst = sum(~isnan(lickBurstStart(:,ind_omit)));
            nnan = sum(isnan(lickBurstStart(:,ind_omit)));
            ind_early_omit = sortlick_ind(1:floor(nburst/4));
            ind_late_omit = sortlick_ind(nburst-floor(nburst/4)+1:end-nnan);
            early_omit_time = mean((lickBurstStart(:,ind_omit(ind_early_omit))-prewin_frames).*(1000./frameRateHz),2);
            late_omit_time = mean((lickBurstStart(:,ind_omit(ind_late_omit))-prewin_frames).*(1000./frameRateHz),2);
        else
            ind_early_omit = [];
            ind_late_omit = [];
            early_omit_time = NaN;
            late_omit_time = NaN;
        end
        if sum(~isnan(lickBurstStart(:,ind_rew)))>6
        figure;
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_early_rew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_early_rew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC), 'k');
        hold on
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew(ind_late_rew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew(ind_late_rew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b');
        scatter((lickBurstStart(:,ind_rew(ind_early_rew))-prewin_frames).*(1000./frameRateHz),-0.5.*ones(1,length(ind_early_rew)),'xk')
        scatter((lickBurstStart(:,ind_rew(ind_late_rew))-prewin_frames).*(1000./frameRateHz),-0.5.*ones(1,length(ind_late_rew)),'xb')
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([-1 inf])
        title([mouse ' ' date '- lick bursts: early (n = ' num2str(length(ind_early_rew)) '); late (n = ' num2str(length(ind_late_rew)) ')'])
        savefig(fullfile(lg_out,img_fn, [img_fn '_cueAlignSpiking_byLickTime.fig']))
        savefig(fullfile(lg_out,img_fn, [img_fn '_cueAlignSpiking_byLickTime.fig']))
        end
        
        pct_precue_burst = length(find((lickBurstStart-prewin_frames).*(1000./frameRateHz)<600))./size(lickBurstStart,2);
        tl = (1-postLick_frames:postLick_frames).*(1000./frameRateHz);
        
        figure; 
        subplot(2,3,1)
        shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents(:,:,ind_rew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
        hold on
        if sum(~isnan(postRew_lickBurstStart(ind_omit)))>4
            shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents(:,:,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents(:,:,ind_omit),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
        end
        xlabel('Time from lick')
        ylabel('Spike rate (Hz)')
        ylim([0 inf])
        title(['Reward (black- ' num2str(sum(~isnan(squeeze(postRew_lickAlignEvents(1,1,ind_rew))))) '); omit (red- ' num2str(sum(~isnan(squeeze(postRew_lickAlignEvents(1,1,ind_omit))))) ')'])
        subplot(2,3,2)
        [sortlick sortlick_ind] = sort(postRew_lickBurstStart(:,ind_rew),'ascend');
        nburst = sum(~isnan(postRew_lickBurstStart(:,ind_rew)));
        nnan = sum(isnan(postRew_lickBurstStart(:,ind_rew)));
        ind_prerew_early_rew = sortlick_ind(1:floor(nburst/4));
        ind_prerew_late_rew = sortlick_ind(nburst-floor(nburst/4)+1:end-nnan);
        early_rew_time = nanmean((postRew_lickBurstStart(:,ind_rew(ind_prerew_early_rew))-prewin_frames-rewDelay_frames).*(1000./frameRateHz),2);
        late_rew_time = nanmean((postRew_lickBurstStart(:,ind_rew(ind_prerew_late_rew))-prewin_frames-rewDelay_frames).*(1000./frameRateHz),2);
        shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents(:,:,ind_rew(ind_prerew_early_rew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents(:,:,ind_rew(ind_prerew_early_rew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
        hold on
        shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents(:,:,ind_rew(ind_prerew_late_rew)),3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents(:,:,ind_rew(ind_prerew_late_rew)),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'b');
        xlabel('Time from lick')
        ylabel('Spike rate (Hz)')
        ylim([0 inf])
        title(['Black- avg = ' num2str(chop(early_rew_time,3)) 'ms; Blue- avg = ' num2str(chop(late_rew_time,3)) 'ms)'])
        subplot(2,3,4)
        shadedErrorBar(tl, nanmean(postRew_lickAlign(:,ind_rew),2).*(1000./frameRateHz), (nanstd(postRew_lickAlign(:,ind_rew),[],2).*(1000./frameRateHz))./sqrt(length(ind_rew)),'k');
        hold on
        if sum(~isnan(postRew_lickBurstStart(ind_omit)))>4
            shadedErrorBar(tl, nanmean(postRew_lickAlign(:,ind_omit),2).*(1000./frameRateHz), (nanstd(postRew_lickAlign(:,ind_omit),[],2).*(1000./frameRateHz))./sqrt(length(ind_omit)),'r');
        end
        xlabel('Time from lick')
        ylabel('Lick rate (Hz)')
        ylim([0 inf])
        subplot(2,3,5)
        shadedErrorBar(tl, nanmean(postRew_lickAlign(:,ind_rew(ind_prerew_early_rew)),2).*(1000./frameRateHz), (nanstd(postRew_lickAlign(:,ind_rew(ind_prerew_early_rew)),[],2).*(1000./frameRateHz))./sqrt(length(ind_rew(ind_prerew_early_rew))),'k');
        hold on
        shadedErrorBar(tl, nanmean(postRew_lickAlign(:,ind_rew(ind_prerew_late_rew)),2).*(1000./frameRateHz), (nanstd(postRew_lickAlign(:,ind_rew(ind_prerew_late_rew)),[],2).*(1000./frameRateHz))./sqrt(length(ind_rew(ind_prerew_late_rew))),'b');
        xlabel('Time from lick')
        ylabel('Lick rate (Hz)')
        ylim([0 inf])
        if id == 3
            subplot(2,3,3)
            shadedErrorBar(tl, nanmean(nanmean(postRew_lickAlignEvents(:,:,ind_unexp),3),2).*(1000./frameRateHz), (nanstd(nanmean(postRew_lickAlignEvents(:,:,ind_unexp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'g');
            xlabel('Time from lick')
            ylabel('Spike rate (Hz)')
            ylim([0 inf])
            title(['Unexpected reward (n = ' num2str(sum(~isnan(squeeze(postRew_lickAlignEvents(1,1,ind_unexp))))) ')'])
            subplot(2,3,6)
            shadedErrorBar(tl, nanmean(postRew_lickAlign(:,ind_unexp),2).*(1000./frameRateHz), (nanstd(postRew_lickAlign(:,ind_unexp),[],2).*(1000./frameRateHz))./sqrt(length(ind_unexp)),'g');
            xlabel('Time from lick')
            ylabel('Lick rate (Hz)')
            ylim([0 inf])
        end
        suptitle([mouse ' ' date '- post reward lick burst aligned spiking'])
        savefig(fullfile(lg_out,img_fn, [img_fn '_postRew_lickAlignSpiking.fig']))
        
        figure; 
        for i = 1:4
            plot(tt,cumsum(sum(lickCueAlign(:,1+((i-1).*floor(nTrials/4)):floor(nTrials/4)+((i-1).*floor(nTrials/4))),2)));
            hold on 
        end
        title([mouse ' ' date '- cumulative licking by quarter session'])
        xlabel('Time from cue')
        ylabel('Cumulative Licks')
        savefig(fullfile(lg_out,img_fn, [img_fn '_cumulativeLicking.fig']))   
        
        tl_rew = (1-postLick_frames:(postLick_frames*2)).*(1000./frameRateHz);
        save(fullfile(lg_out,img_fn, [img_fn '_cueAlignLick.mat']), 'firstPostRewLickFrame', 'lastPreRewLickFrame', 'tl_rew', 'firstPostRew_lickAlignEvents', 'lastPreRew_lickAlignEvents', 'lickCueAlign', 'lickBurstStart', 'lickCounterVals', 'lickSearch_frames', 'lickDelay_frames', 'lickTC', 'postwin_frames', 'prewin_frames', 'frameRateHz', 'tt', 'ind_early_rew', 'ind_late_rew', 'ind_early_omit', 'ind_late_omit', 'early_rew_time', 'late_rew_time', 'early_omit_time', 'late_omit_time','pct_precue_burst', 'postRew_lickAlignEvents', 'postLick_frames', 'postRew_lickBurstStart','tl', 'ind_prerew_early_rew', 'ind_prerew_late_rew','postRew_lickAlign')

        figure;
        subplot(2,2,1)
        shadedErrorBar(tl_rew, nanmean(nanmean(lastPreRew_lickAlignEvents(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(lastPreRew_lickAlignEvents(:,:,ind_rew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
        hold on
        if length(ind_omit)>5
            shadedErrorBar(tl_rew, nanmean(nanmean(lastPreRew_lickAlignEvents(:,:,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(lastPreRew_lickAlignEvents(:,:,ind_omit),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
        end
        if length(ind_unexp)>5
            shadedErrorBar(tl_rew, nanmean(nanmean(lastPreRew_lickAlignEvents(:,:,ind_unexp),3),2).*(1000./frameRateHz), (nanstd(nanmean(lastPreRew_lickAlignEvents(:,:,ind_unexp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'g');
        end
        title('Last lick before reward')
        xlabel('Time from lick (ms)')
        ylabel('Spike rate (Hz)')
        ylim([0 6])
        subplot(2,2,3)
        shadedErrorBar(tl_rew, nanmean(lastPreRew_lickAlign(:,ind_rew),2).*(1000./frameRateHz), (nanstd(lastPreRew_lickAlign(:,ind_rew),[],2).*(1000./frameRateHz))./sqrt(length(ind_rew)),'k');
        hold on
        if length(ind_omit)>5
            shadedErrorBar(tl_rew, nanmean(lastPreRew_lickAlign(:,ind_omit),2).*(1000./frameRateHz), (nanstd(lastPreRew_lickAlign(:,ind_omit),[],2).*(1000./frameRateHz))./sqrt(length(ind_omit)),'r');
            omit_str = [' Omit- ' num2str(length(intersect(ind_omit,preRewTrials))) ' trials; '];
        else
            omit_str = [];
        end
        if length(ind_unexp)>5
            shadedErrorBar(tl_rew, nanmean(lastPreRew_lickAlign(:,ind_unexp),2).*(1000./frameRateHz), (nanstd(lastPreRew_lickAlign(:,ind_unexp),[],2).*(1000./frameRateHz))./sqrt(length(ind_unexp)),'g');        
            unexp_str = ['Unexp- ' num2str(length(intersect(ind_unexp,preRewTrials))) ' trials'];
        else
            unexp_str = [];
        end
        xlabel('Time from lick (ms)')
        ylabel('Lick rate (Hz)')
        ylim([0 30])
        title(['Rew- ' num2str(length(intersect(ind_rew,preRewTrials))) ' trials; ' omit_str unexp_str])        
        subplot(2,2,2)
        shadedErrorBar(tl_rew, nanmean(nanmean(firstPostRew_lickAlignEvents(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(firstPostRew_lickAlignEvents(:,:,ind_rew),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
        hold on
        if length(ind_omit)>5
            shadedErrorBar(tl_rew, nanmean(nanmean(firstPostRew_lickAlignEvents(:,:,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(firstPostRew_lickAlignEvents(:,:,ind_omit),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
        end
        if length(ind_unexp)>5
            shadedErrorBar(tl_rew, nanmean(nanmean(firstPostRew_lickAlignEvents(:,:,ind_unexp),3),2).*(1000./frameRateHz), (nanstd(nanmean(firstPostRew_lickAlignEvents(:,:,ind_unexp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'g');
        end   
        title('First lick after reward')
        xlabel('Time from lick (ms)')
        ylabel('Spike rate (Hz)')
        ylim([0 6])
        subplot(2,2,4)
        shadedErrorBar(tl_rew, nanmean(firstPostRew_lickAlign(:,ind_rew),2).*(1000./frameRateHz), (nanstd(firstPostRew_lickAlign(:,ind_rew),[],2).*(1000./frameRateHz))./sqrt(length(ind_rew)),'k');
        hold on
        if length(ind_omit)>5
            shadedErrorBar(tl_rew, nanmean(firstPostRew_lickAlign(:,ind_omit),2).*(1000./frameRateHz), (nanstd(firstPostRew_lickAlign(:,ind_omit),[],2).*(1000./frameRateHz))./sqrt(length(ind_omit)),'r');
            omit_str = [' Omit- ' num2str(length(intersect(ind_omit,postRewTrials))) ' trials; '];
        else
            omit_str = [];
        end
        if length(ind_unexp)>5
            shadedErrorBar(tl_rew, nanmean(firstPostRew_lickAlign(:,ind_unexp),2).*(1000./frameRateHz), (nanstd(firstPostRew_lickAlign(:,ind_unexp),[],2).*(1000./frameRateHz))./sqrt(length(ind_unexp)),'g');        
            unexp_str = ['Unexp- ' num2str(length(intersect(ind_unexp,postRewTrials))) ' trials'];
        else
            unexp_str = [];
        end
        xlabel('Time from lick (ms)')
        ylabel('Lick rate (Hz)')
        ylim([0 30])
        title(['Rew- ' num2str(length(intersect(ind_rew,postRewTrials))) ' trials; ' omit_str unexp_str])
        suptitle([mouse ' ' date '- Licks relative to reward'])
        savefig(fullfile(lg_out,img_fn, [img_fn '_lastVsFirstLick.fig'])) 
    end
    
    close all
end
    