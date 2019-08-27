clear all
plotDetect = 0;
negAdjust = 0;
for id = [1,2]%[1,2,3,5,6]
    %close all
    %CRP_expt_list_Crus_jh;
    CRP_expt_list_LS_jh;
    jake_out = 'Y:\home\jake\Analysis\Cue_reward_pairing_analysis\CC_analysis_2P';
    behav_dir = 'Y:\home\jake\Data\WidefieldImaging\GCaMP\behavior';
    nexp = size(expt(id).date,1);
    
%     if id==1  %Crus
%         exp_subset= [5:7]; % Day1
%     elseif id==2
%         exp_subset = [6:8]; % PL OR
%     elseif id==3
%         exp_subset = [4:6]; % UR
%     elseif id==4
%         exp_subset = [1:3]; %  OT OR
%     elseif id==5
%         exp_subset = [1:3]; % OT UR
%     end
     
    if id==1  %LS
        exp_subset=  [18:20]; %
    elseif id==2
        exp_subset = [19:20]; %
    elseif id==3
        exp_subset = [7:9]; %
    elseif id==5
        exp_subset = [3:5]; %
    elseif id==6
        exp_subset = [1:2]; %
    end    
    
    for iexp = exp_subset  %1:nexp
        mouse = strtrim(expt(id).mouse(iexp,:));
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];
        load(fullfile(jake_out,img_fn, [img_fn '_ROI_TCs.mat']))
        nIC = size(tc_avg,2);

        opt.thresh = 1.7;
        opt.normalization = 1;
        opt.deconvtau = 0;
        opt.dt = 1;

%         if str2num(mouse)>800
%             mouse_name = mouse(2:3);
%         else
%             mouse_name = mouse;
%         end
%         
        if str2num(mouse)>800 & str2num(mouse) < 1000
            mouse_name = mouse(2:3);
        elseif str2num(mouse)<100 | str2num(mouse) > 1000
            mouse_name = mouse;
        end
        %% correct for skipped frames
        fns = dir(fullfile(behav_dir, ['data-i' mouse '-' date '*']));
        if length(fns)>1
            nfns{id}(iexp) = length(fns);
            fn(id).tr{iexp}= [];
            for i = 1:length(fns)
                load(fullfile(behav_dir,fns(i).name));
                fn(id).tr{iexp} = [fn(id).tr{iexp} size(input.gratingContrast,2)];
            end
            fns = fns(end); %assumes the final behavioral file collected is the correct one
        else
            nfns{id}(iexp) = 1;  
        end
        load(fullfile(behav_dir,fns.name));  %laods bx data
        nTrials = length(input.cTargetOn);
        input.counterValues_orig = input.counterValues; 
        input.counterTimesUs_orig = input.counterTimesUs;
        input.cTargetOn_orig = input.cTargetOn;
        val_reset = 0;
        for it = 1:nTrials
            diff_val = diff(input.counterValues{it});
            if length(unique(diff_val)) == 1
                continue
            elseif length(unique(diff_val)) == 2 & sum(ismember(unique(diff_val),2)) == 1
                val_reset = 1;
                ind = find(diff_val== 2);
                for i = 1:length(ind)
                    d = input.counterTimesUs{it}(ind(i))-input.counterTimesUs{it}(ind-1);
                    if d>50000 & d<70000
                        val_reset = 1;
                        input.counterValues{it} = [input.counterValues{it}(1:ind(i)-1) input.counterValues{it}(1:ind(i)-1)+1 input.counterValues{it}(ind(i):end)];
                        input.counterValues{it} = [input.counterTimesUs{it}(1:ind(i)-1) input.counterTimesUs{it}(1:ind(i)-1)+33000 input.counterTimesUs{it}(ind(i):end)];
                        ind = ind+1;
                    end
                end
            elseif length(unique(diff_val)) > 2
                val_reset = 1;
                ind = find(diff(input.counterTimesUs{it}./1000)<25);
                ind = ind(1)-2:ind(end)+2;
                timeDiff = input.counterTimesUs{it}(ind(end))-input.counterTimesUs{it}(ind(1));
                frameDiff = input.counterValues{it}(ind(end))-input.counterValues{it}(ind(1));
                nf = round((timeDiff./1000)*frameRateHz/1000);
                input.counterValues{it} = [input.counterValues{it}(1:ind(1)-1) bsxfun(@plus,input.counterValues{it}(ind(1)-1), 1:nf) [input.counterValues{it}(ind(end):end)-frameDiff+nf]];
                input.counterTimesUs{it} = [input.counterTimesUs{it}(1:ind(1)-1) [repmat(input.counterTimesUs{it}(ind(1)-1),[1 nf]) + (33000:33000:33000*nf)] input.counterTimesUs{it}(ind(end):end)];
                if input.cTargetOn{it}>ind(end)
                    input.cTargetOn{it} = input.cTargetOn{it}-frameDiff+nf;
                elseif input.cTargetOn{it}>ind(1)
                    input.cTargetOn{it} = NaN;
                end
                for it2 = it+1:nTrials
                    input.counterValues{it2} = bsxfun(@plus, input.counterValues_orig{it2},-frameDiff+nf);
                    input.cTargetOn{it2} = bsxfun(@plus, input.cTargetOn_orig{it2}, -frameDiff+nf);
                end
                diff_val = diff(input.counterValues{it});
                if length(unique(diff_val)) == 1
                    continue
                elseif length(unique(diff_val)) == 2 & sum(ismember(unique(diff_val),2))
                    ind = find(diff_val== 2);
                    for i = 1:length(ind)
                        d = input.counterTimesUs{it}(ind(i))-input.counterTimesUs{it}(ind-1);
                        if d>50000 & d<70000
                            val_reset = 1;
                            input.counterValues{it} = [input.counterValues{it}(1:ind(i)-1) input.counterValues{it}(1:ind(i)-1)+1 input.counterValues{it}(ind(i):end)];
                            input.counterValues{it} = [input.counterTimesUs{it}(1:ind(i)-1) input.counterTimesUs{it}(1:ind(i)-1)+33000 input.counterTimesUs{it}(ind(i):end)];
                            ind = ind+1;
                        end
                    end
                end
            end
        end
        if val_reset == 0
            rmfield(input,{'counterValues_orig','counterTimesUs_orig','cTargetOn'});
            fprintf(' counter not corrected \n')
        else
            fprintf(' corrected counter \n')
        end

        %% detect events

        tc_avg_temp = tc_avg;
        cTargetOn = celleqel2mat_padded(input.cTargetOn);
        if expt(id).ttl(iexp)
            cd(fullfile('Y:\home\jake\Data\2P_imaging', [date '_img' mouse_name],['img' mouse_name]))
            load(['img' mouse_name '_000_' run '_realtime.mat'])
            %if ttl_log is faulty then identify pockel transitions via changes in F
            if length(unique(ttl_log)) == 1 | sum(ttl_log) < 0.05*(length(ttl_log))
                clear ttl_log
                ttl_log = find_pockel_tc(tc_avg, 1)';
            end
            ttl_log = ttl_log(1:size(tc_avg,1),:);
            ttl_trans = find(abs(diff(ttl_log))); 
            n = length(ttl_trans);
            for i = 1:n
                ttl_log(ttl_trans(i)-4:ttl_trans(i)+4,:) = 0;
            end
            ind_ttl = intersect(1:size(tc_avg,1), find(ttl_log == 0));
            tc_avg(ind_ttl,:) = NaN;
            ind_tc = find(~isnan(tc_avg(:,1)));
            tc_avg_temp(find(isnan(tc_avg(:,1))),:) = [];
        end
        range = 5501:5900;
        
        tc_diff = diff(tc_avg_temp,[],1);
        all_events = zeros(size(tc_avg_temp));
        if plotDetect
            figure;
            start = 1;
            IC_list = [];
        end
        for ic = 1:nIC
            [~, iti_ind, ~] = CellsortFindspikes(tc_diff(:,ic), opt.thresh, opt.dt, opt.deconvtau, opt.normalization);
            iti_ind(find(diff(iti_ind)==1)+1)=[];
            all_events(iti_ind+1,ic) = 1;
            if negAdjust
                ind_events = find(all_events(:,ic));
                for i = 1:length(ind_events)
                    if ind_events(i)< 50
                        if tc_avg_temp(ind_events(i)-1,ic) < mean(tc_avg_temp(1:ind_events(i)+50,ic),1) - 2.*std(tc_avg_temp(1:ind_events(i)+50,ic),[],1)
                            all_events(ind_events(i),ic) = 0;
                        end
                    elseif size(tc_avg_temp,1)-ind_events(i)<50
                        if tc_avg_temp(ind_events(i)-1,ic) < mean(tc_avg_temp(ind_events(i)-50:end,ic),1) - 2.*std(tc_avg_temp(ind_events(i)-50:end,ic),[],1)
                            all_events(ind_events(i),ic) = 0;
                        end
                    else
                        if tc_avg_temp(ind_events(i)-1,ic) < mean(tc_avg_temp(ind_events(i)-50:ind_events(i)+50,ic),1) - 2.*std(tc_avg_temp(ind_events(i)-50:ind_events(i)+50,ic),[],1)
                            all_events(ind_events(i),ic) = 0;
                        end
                    end
                end
            end
            if plotDetect
                if start>9
                    savefig(fullfile(jake_out,img_fn, [img_fn '_Cell#' num2str(IC_list(1)) '-' num2str(IC_list(end)) '_diff_findSpikes.fig']))
                    figure; start = 1;
                    IC_list = [];
                end

                subplot(3,3,start)

                plot((tc_avg_temp(range,ic)-nanmean(tc_avg_temp(range,ic),1))./max(tc_avg_temp(range,ic),[],1)); 
                hold on; plot(all_events(range,ic));
                title(['cell#' num2str(ic)])
                start = start+1;
                IC_list = [IC_list ic];
            end
        end
        if plotDetect
            savefig(fullfile(jake_out,img_fn, [img_fn '_Cell#' num2str(IC_list(1)) '-' num2str(IC_list(end)) '_diff_findSpikes.fig']))
        end
        if expt(id).ttl(iexp)
            all_events_temp = nan(size(tc_avg));
            all_events_temp(ind_tc,:) = all_events;
            all_events = all_events_temp;
        end

        save(fullfile(jake_out,img_fn, [img_fn '_ROI_spikes.mat']), 'all_events', 'opt', 'tc_diff')

        if expt(id).ttl(iexp)
            nFrames = length(ind_tc);
        else
            nFrames = size(all_events,1);
        end
        frameRateHz = double(input.frameRateHz);
        spikeRate = (nansum(all_events,1)./nFrames).*frameRateHz; 
        figure; hist(spikeRate)
        xlim([0 5]); xlabel('Spike rate'); ylabel('Cells'); 
        title([date ' ' mouse ' spike rates'])
        savefig(fullfile(jake_out,img_fn, [img_fn '_spikeRate.fig']))

        
        %% align events    
        
        omitRewTrial = celleqel2mat_padded(input.tRewardOmissionTrial);
        unexpRewTrial = celleqel2mat_padded(input.tDoNoStimulusChange);
        nTrials = length(cTargetOn);
        prewin_frames = round(1500./frameRateHz);
        postwin_frames = round(3000./frameRateHz);
        targetAlign_tc = nan(prewin_frames+postwin_frames,nIC,nTrials);
        targetAlign_events = nan(prewin_frames+postwin_frames,nIC,nTrials);
        nFrames = size(tc_avg,1);
        
        for itrial = 1:nTrials
            if cTargetOn(itrial)+postwin_frames-1 <= nFrames %& input.counterValues{itrial}(end)-cTargetOn(itrial) > postwin_frames
                targetAlign_tc(:,:,itrial) = tc_avg(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
                targetAlign_events(:,:,itrial) = all_events(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:);
                if id == 2 & (iexp == 12 | iexp ==14)
                    targetAlign_tc(prewin_frames+3:prewin_frames+5,:,itrial) = NaN;
                    targetAlign_events(prewin_frames+3:prewin_frames+5,:,itrial) = NaN;
                end
            end
        end

        targetAlignF = nanmean(targetAlign_tc(1:prewin_frames,:,:),1);
        targetAligndFoverF = (targetAlign_tc-repmat(targetAlignF, [size(targetAlign_tc,1),1,1])./repmat(targetAlignF, [size(targetAlign_tc,1),1,1]));

        ind_omit = find(omitRewTrial);
        ind_unexp = find(unexpRewTrial);
        ind_rew = intersect(find(omitRewTrial == 0),find(unexpRewTrial == 0));
        if intersect(ind_omit, ind_unexp)
            x = ismember(ind_omit,intersect(ind_omit, ind_unexp));
            ind_omit(x) = [];
            x = ismember(ind_unexp,intersect(ind_omit, ind_unexp));
            ind_unexp(x) = [];
        end
        tt = (-prewin_frames:postwin_frames-1).*(1000./frameRateHz);
        rewardDelayDurationMs = double(max(celleqel2mat_padded(input.tRewardDelayDurationMs),[],2));
        reactTimesMs = double(mean(celleqel2mat_padded(input.reactTimesMs),2));
        delayTimeMs = reactTimesMs+rewardDelayDurationMs;


        figure;
        subplot(3,1,1)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_rew),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_rew),3),[],2)./sqrt(nIC));
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        title('Reward')
        subplot(3,1,2)
        if length(ind_omit>5)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_omit),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_omit),3),[],2)./sqrt(nIC),'r');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        title('Omit')
        end
        subplot(3,1,3)
        if length(ind_unexp>5)
        shadedErrorBar(tt, nanmean(nanmean(targetAligndFoverF(:,:,ind_unexp),3),2), nanstd(nanmean(targetAligndFoverF(:,:,ind_unexp),3),[],2)./sqrt(nIC),'g');
        xlabel('Time from cue')
        ylabel('Avg dF/F')
        title('Unexpected reward')
        end
        %suptitle([date ' ' mouse])
        savefig(fullfile(jake_out,img_fn, [img_fn '_cueAlign_dFoverF.fig']))

        figure;
        subplot(3,1,1)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew),3),[],2)./sqrt(nIC)).*(1000./frameRateHz));
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Reward')
        if length(ind_omit>5)
        subplot(3,1,2)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'r');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omit')
        end
        if length(ind_unexp>5)
        subplot(3,1,3)
        shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_unexp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_unexp),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'g');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Unexpected reward')
        end    
        %suptitle([date ' ' mouse])
        savefig(fullfile(jake_out,img_fn, [img_fn '_cueAlign_events_Hz.fig']))
        
        if length(ind_omit>10)
            figure;
            subplot(2,1,1)
            thresh = mean(diff(ind_omit));
            ind_omit_short = find(diff(ind_omit)<thresh)+1;
            ind_omit_long = find(diff(ind_omit)>thresh)+1;
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_short)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_short)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_long)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit(ind_omit_long)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Short (black n = ' num2str(length(ind_omit_short)) ') vs long (blue n = ' num2str(length(ind_omit_long)) ') interval between omits'])
            ind_rew_preomit = ind_omit-1;
            if find(ind_rew_preomit==0)
              ind_rew_preomit(find(ind_rew_preomit==0)) = [];
            end
            ind_rew_postomit = ind_omit+1;
            if find(ind_rew_postomit)
              ind_rew_postomit(find(ind_rew_postomit>nTrials)) = [];
            end
            subplot(2,1,2)
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew_preomit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_preomit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_rew_postomit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_postomit),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Reward trial pre (black n = ' num2str(length(ind_rew_preomit)) ') vs post (blue n = ' num2str(length(ind_rew_postomit)) ') omit trials'])
            %suptitle([date ' ' mouse])
            savefig(fullfile(jake_out,img_fn, [img_fn '_omitByInterval.fig']))
        else
            ind_omit_short = [];
            ind_omit_long = [];
            ind_rew_preomit = [];
            ind_rew_postomit = [];
        end
        if length(ind_unexp>10)
            figure;
            thresh = mean(diff(ind_unexp));
            ind_unexp_short = find(diff(ind_unexp)<thresh)+1;
            ind_unexp_long = find(diff(ind_unexp)>thresh)+1;
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_short)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_short)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_long)),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_unexp(ind_unexp_long)),3),[],2)./sqrt(nIC)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title([date ' ' mouse '- Short (black n = ' num2str(length(ind_unexp_short)) ') vs long (blue n = ' num2str(length(ind_unexp_long)) ') interval between unexpected reward'])
            savefig(fullfile(jake_out,img_fn, [img_fn '_unexpByInterval.fig']))
        else
            ind_unexp_short = [];
            ind_unexp_long = [];
        end
        
        figure;
        n = floor(nTrials./3);
        start = 1;
        for i = 1:3
            subplot(3,2,start)
            ind_rew_temp = intersect(ind_rew,1+(i-1).*n:i*n);
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_rew_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_rew_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'k');
            hold on
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
            vline(600)
            if length(ind_omit)>=10
                subplot(3,2,start+1)
                ind_omit_temp = intersect(ind_omit,1+(i-1).*n:i*n);
                shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,ind_omit_temp),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,ind_omit_temp),3),[],2).*(1000./frameRateHz))./sqrt(nIC),'r');
                ylim([0 5])
                xlabel('Time from cue')
                ylabel('Spike rate (Hz)')
                title(['Trials ' num2str(1+(i-1).*n) ':' num2str(i*n)])
                vline(600)
            end
            start = start+2;
        end
        %suptitle([date ' ' mouse '- Reward (black), Omit (red)'])
        savefig(fullfile(jake_out,img_fn, [img_fn '_repsByTrial.fig']))
        
        if unique(expt(id).name == 'Crus')
            load(fullfile(jake_out,img_fn, [img_fn '_splitImage.mat']))
            figure;
            indL = find(maskCat==1);
            indR = find(maskCat==2);
            subplot(2,2,1)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indL,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indL,ind_rew),3),[],2).*(1000./frameRateHz))./sqrt(length(indL)),'k');
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Reward- Left side- n=' num2str(length(indL))])
            vline(600)
            subplot(2,2,2)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indR,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indR,ind_rew),3),[],2).*(1000./frameRateHz))./sqrt(length(indR)),'k');
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            vline(600)
            title(['Reward- Right side- n=' num2str(length(indR))])
            subplot(2,2,3)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indL,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indL,ind_omit),3),[],2).*(1000./frameRateHz))./sqrt(length(indL)),'r');
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['Omit- Left side- n=' num2str(length(indL))])
            vline(600)
            subplot(2,2,4)
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indR,ind_omit),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indR,ind_omit),3),[],2).*(1000./frameRateHz))./sqrt(length(indR)),'r');
            ylim([0 5])
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['omit- Right side- n=' num2str(length(indR))])
            vline(600)
            %suptitle([date ' ' mouse '- Reward (black), Omit (red)'])
            savefig(fullfile(jake_out,img_fn, [img_fn '_repsByCrus.fig']))
        end

        
        save(fullfile(jake_out,img_fn, [img_fn '_targetAlign.mat']), 'ind_rew', 'ind_omit', 'ind_unexp', 'targetAlign_events', 'targetAligndFoverF', 'prewin_frames', 'postwin_frames', 'tt', 'frameRateHz', 'ind_omit_short','ind_omit_long','ind_unexp_short','ind_unexp_long','ind_rew_preomit','ind_rew_postomit')
        save(fullfile(jake_out,img_fn, [img_fn '_input.mat']), 'input')
        %close all
    end
end
