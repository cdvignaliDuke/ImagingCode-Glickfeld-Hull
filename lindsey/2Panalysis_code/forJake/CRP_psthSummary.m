clear all
close all
CRP_expt_list_Crus
plotAvg = 0;
doSpr2017sep = 0;
area_list = {'C1','C2','LS'};
for id = 3
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
    if doSpr2017sep
        all_rew_6 = [];
        all_omit_6 = [];
        all_unexp_6 = [];
        all_early_rew_6 = [];
        all_late_rew_6 = [];
        all_early_omit_6 = [];
        all_late_omit_6 = [];
        all_rew_6_df = [];
        all_omit_6_df = [];
        all_unexp_6_df = [];
        all_lick_rew_6 = [];
        all_lick_omit_6 = [];
        all_lick_unexp_6 = [];
        all_early_rew_time_6 = zeros(1,nexp);
        all_late_rew_time_6 = zeros(1,nexp);
        all_early_omit_time_6 = zeros(1,nexp);
        all_late_omit_time_6 = zeros(1,nexp);
        all_short_omit_6 = [];
        all_long_omit_6 = [];
        all_short_unexp_6 = [];
        all_long_unexp_6 = [];
        all_preomit_6 = [];
        all_postomit_6 = [];
        mouse_str_6 = [];
    end
    all_pct_precue_burst = [];
    mouse_str = [];
    all_area_id = [];
    expt_areas = zeros(length(area_list),nexp);
    for iexp = 1:nexp
        mouse = expt(id).mouse(iexp,:);
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];
        mouse_str = [mouse_str '_' mouse];
        
        load(fullfile(lg_out,img_fn, [img_fn '_targetAlign.mat']))
        load(fullfile(lg_out,img_fn, [img_fn '_cueAlignLick.mat']))
        if exist(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']))
            load(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']))
            area_temp = zeros(1,size(targetAligndFoverF,2));
            for i = 1:2
                x = find(strcmp(area_list, expt(id).areas{iexp}(i)));
                ind = find(maskCat==i);
                area_temp(1,ind) = x;
                expt_areas(x,iexp) = 1;
            end
            all_area_id = [all_area_id area_temp];
        end
        
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
            postrew_lick_omit_avg = [];
            postrew_lick_omit_lick_avg = [];
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

        if ~exist(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id)])
            mkdir(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id)])
        end
        
        if plotAvg
        figure;
        IC_list = [];
        start = 1;
        for ic = 1:nIC
            if start > 49
                suptitle([date ' ' mouse '- Cell#' num2str(IC_list(1)) ':' num2str(IC_list(end)) '; Trials: ' num2str(n_rew) ' Rew, ' num2str(n_omit) ' Omit, '  num2str(n_unexp) ' Unexp']) 
                print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\' img_fn '_Cell#' num2str(IC_list(1)) '-' num2str(IC_list(end)) '_psth.pdf'],'-dpdf', '-fillpage')
                figure;
                start = 1;
                IC_list = [];
            end
            subplot(7,7,start)
            hold on
            if length(ind_omit>1)
                shadedErrorBar(tt, omit_avg(:,ic).*(1000./frameRateHz), (nanstd(targetAlign_events(:,ic,ind_omit),[],3)./sqrt(length(ind_omit))).*(1000./frameRateHz),'r');
                hold on
            end
            if length(ind_unexp>1)
                shadedErrorBar(tt, unexp_avg(:,ic).*(1000./frameRateHz), (nanstd(targetAlign_events(:,ic,ind_unexp),[],3)./sqrt(length(ind_unexp))).*(1000./frameRateHz),'g');
                hold on
                vline(1100, 'b')
            end
            shadedErrorBar(tt, rew_avg(:,ic).*(1000./frameRateHz), (nanstd(targetAlign_events(:,ic,ind_rew),[],3)./sqrt(length(ind_rew))).*(1000./frameRateHz),'k');
            xlim([-500 2000])
            vline([0 600], 'b')
            start = start+1;
            IC_list = [IC_list ic];
        end
        suptitle([date ' ' mouse '- Cell#' num2str(IC_list(1)) ':' num2str(IC_list(end)) '; Trials: ' num2str(n_rew) ' Rew, ' num2str(n_omit) ' Omit, '  num2str(n_unexp) ' Unexp']) 
        print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\' img_fn '_Cell#' num2str(IC_list(1)) '-' num2str(IC_list(end)) '_psth.pdf'],'-dpdf', '-fillpage')
        
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
        all_pct_precue_burst = [all_pct_precue_burst pct_precue_burst];
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

        if iexp <= 6 & doSpr2017sep
            all_rew_6 = [all_rew_6 nanmean(targetAlign_events(:,:,ind_rew),3)];
            all_omit_6 = [all_omit_6 nanmean(targetAlign_events(:,:,ind_omit),3)];
            all_unexp_6 = [all_unexp_6 nanmean(targetAlign_events(:,:,ind_unexp),3)];
            all_early_rew_6 = [all_early_rew_6 nanmean(targetAlign_events(:,:,ind_rew(ind_early_rew)),3)];
            all_late_rew_6 = [all_late_rew_6 nanmean(targetAlign_events(:,:,ind_rew(ind_late_rew)),3)];
            all_early_omit_6 = [all_early_omit_6 nanmean(targetAlign_events(:,:,ind_rew(ind_early_omit)),3)];
            all_late_omit_6 = [all_late_omit_6 nanmean(targetAlign_events(:,:,ind_rew(ind_late_omit)),3)];
            all_short_omit_6 = [all_short_omit short_omit_avg];
            all_long_omit_6 = [all_long_omit_6 long_omit_avg];
            all_short_unexp_6 = [all_short_unexp_6 short_unexp_avg];
            all_long_unexp_6 = [all_long_unexp_6 long_unexp_avg];
            all_preomit_6 = [all_preomit_6 preomit_avg];
            all_postomit_6 = [all_postomit_6 postomit_avg];
            all_rew_6_df = [all_rew_6_df nanmean(targetAligndFoverF(:,:,ind_rew),3)];
            all_omit_6_df = [all_omit_6_df nanmean(targetAligndFoverF(:,:,ind_omit),3)];
            all_unexp_6_df = [all_unexp_6_df nanmean(targetAligndFoverF(:,:,ind_unexp),3)];
            all_lick_rew_6 = [all_lick_rew_6 nanmean(lickCueAlign(:,ind_rew),2)];
            all_lick_omit_6 = [all_lick_omit_6 nanmean(lickCueAlign(:,ind_omit),2)];
            all_lick_unexp_6 = [all_lick_unexp_6 nanmean(lickCueAlign(:,ind_unexp),2)];
            all_early_rew_time_6 = [all_early_rew_time_6 early_rew_time];
            all_late_rew_time_6 = [all_late_rew_time_6 late_rew_time];
            all_early_omit_time_6 = [all_early_omit_time_6 early_omit_time];
            all_late_omit_time_6 = [all_late_omit_time_6 late_omit_time];
            mouse_str_6 = [mouse_str_6 '_' mouse];
        end
    end

    totIC = size(all_rew,2);
    figure;
    subplot(3,3,1)
    shadedErrorBar(tt, nanmean(all_rew_df,2), (nanstd(all_rew_df,[],2)./sqrt(totIC)),'k');
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
    shadedErrorBar(tt, nanmean(all_rew,2).*(1000./frameRateHz), (nanstd(all_rew,[],2)./sqrt(totIC)).*(1000./frameRateHz),'k');
    xlabel('Time from cue')
    ylabel('Spike rate (Hz)')
    title('Reward')
    ylim([0 3])
    xlim([-500 2000])
    vline([600], 'b')
    if id == 4
        vline([1100], 'c')
    end
    subplot(3,3,3)
    shadedErrorBar(tt, nanmean(all_lick_rew,2).*(1000./frameRateHz), (nanstd(all_lick_rew,[],2)./sqrt(nexp)).*(1000./frameRateHz),'k');
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
    shadedErrorBar(tt, nanmean(all_omit_df,2), (nanstd(all_omit_df,[],2)./sqrt(totIC)),'r');
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
    shadedErrorBar(tt, nanmean(all_omit,2).*(1000./frameRateHz), (nanstd(all_omit,[],2)./sqrt(totIC)).*(1000./frameRateHz),'r');
    xlabel('Time from cue')
    ylabel('Spike rate (Hz)')
    title('Omission')
    ylim([0 3])
    xlim([-500 2000])
    vline([600], 'b')
    if id == 4
        vline([1100], 'c')
    end
    subplot(3,3,6)
    shadedErrorBar(tt, nanmean(all_lick_omit,2).*(1000./frameRateHz), (nanstd(all_lick_omit,[],2)./sqrt(nexp)).*(1000./frameRateHz),'r');
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
    shadedErrorBar(tt, nanmean(all_unexp_df,2), (nanstd(all_unexp_df,[],2)./sqrt(totIC)),'g');
    xlabel('Time from cue')
    ylabel('dF/F')
    title('Unexpected Reward')
    ylim([-.05 0.1])
    xlim([-500 2000])
    vline([600], 'b')
    subplot(3,3,8)
    shadedErrorBar(tt, nanmean(all_unexp,2).*(1000./frameRateHz), (nanstd(all_unexp,[],2)./sqrt(totIC)).*(1000./frameRateHz),'g');
    xlabel('Time from cue')
    ylabel('Spike rate (Hz)')
    title('Unexpected Reward')
    ylim([0 3])
    xlim([-500 2000])
    vline([600], 'b')
    subplot(3,3,9)
    shadedErrorBar(tt, nanmean(all_lick_unexp,2).*(1000./frameRateHz), (nanstd(all_lick_unexp,[],2)./sqrt(nexp)).*(1000./frameRateHz),'g');
    xlabel('Time from cue')
    ylabel('Lick rate (Hz)')
    title('Unexpected Reward')
    ylim([0 10])
    xlim([-500 2000])
    vline([600], 'b')
    suptitle(['Day ' num2str(id) ': n= ' num2str(nexp) ' mice, ' num2str(totIC) ' dendrites']);
    savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary' mouse_str '.fig'])    
    
    if size(all_area_id,2)>10
        totExp = sum(expt_areas,2);
        for i = 1:size(area_list,2)
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
    end
            
    if doSpr2017sep
        totIC_6 = size(all_rew_6,2);
        figure;
        subplot(3,3,1)
        shadedErrorBar(tt, nanmean(all_rew_6_df,2), (nanstd(all_rew_6_df,[],2)./sqrt(totIC_6)),'k');
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
        shadedErrorBar(tt, nanmean(all_rew_6,2).*(1000./frameRateHz), (nanstd(all_rew_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'k');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Reward')
        ylim([0 3])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(3,3,3)
        shadedErrorBar(tt, nanmean(all_lick_rew_6,2).*(1000./frameRateHz), (nanstd(all_lick_rew_6,[],2)./sqrt(6)).*(1000./frameRateHz),'k');
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
        shadedErrorBar(tt, nanmean(all_omit_6_df,2), (nanstd(all_omit_6_df,[],2)./sqrt(totIC_6)),'r');
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
        shadedErrorBar(tt, nanmean(all_omit_6,2).*(1000./frameRateHz), (nanstd(all_omit_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'r');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omission')
        ylim([0 3])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        subplot(3,3,6)
        shadedErrorBar(tt, nanmean(all_lick_omit_6,2).*(1000./frameRateHz), (nanstd(all_lick_omit_6,[],2)./sqrt(6)).*(1000./frameRateHz),'r');
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
        shadedErrorBar(tt, nanmean(all_unexp_6_df,2), (nanstd(all_unexp_6_df,[],2)./sqrt(totIC_6)),'g');
        xlabel('Time from cue')
        ylabel('dF/F')
        title('Unexpected Reward')
        ylim([-.05 .1])
        xlim([-500 2000])
        vline([600], 'b')
        subplot(3,3,8)
        shadedErrorBar(tt, nanmean(all_unexp_6,2).*(1000./frameRateHz), (nanstd(all_unexp_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'g');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Unexpected Reward')
        ylim([0 3])
        xlim([-500 2000])
        vline([600], 'b')
        subplot(3,3,9)
        shadedErrorBar(tt, nanmean(all_lick_unexp_6,2).*(1000./frameRateHz), (nanstd(all_lick_unexp_6,[],2)./sqrt(6)).*(1000./frameRateHz),'g');
        xlabel('Time from cue')
        ylabel('Lick rate (Hz)')
        title('Unexpected Reward')
        ylim([0 10])
        xlim([-500 2000])
        vline([600], 'b')
        suptitle(['Day ' num2str(id) ': n= ' num2str(6) ' mice, ' num2str(totIC_6) ' dendrites']);
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_spring2017mice' mouse_str_6 '.fig'])    
    end
    
    figure;
    subplot(2,3,1)
    shadedErrorBar(tt, nanmean(all_earlytrial_rew_df,2), (nanstd(all_earlytrial_rew_df,[],2)./sqrt(totIC)),'k');
    hold on
    shadedErrorBar(tt, nanmean(all_latetrial_rew_df,2), (nanstd(all_latetrial_rew_df,[],2)./sqrt(totIC)),'b');
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
    shadedErrorBar(tt, nanmean(all_earlytrial_rew,2).*(1000./frameRateHz), (nanstd(all_earlytrial_rew,[],2)./sqrt(totIC)).*(1000./frameRateHz),'k');
    hold on
    shadedErrorBar(tt, nanmean(all_latetrial_rew,2).*(1000./frameRateHz), (nanstd(all_latetrial_rew,[],2)./sqrt(totIC)).*(1000./frameRateHz),'b');
    xlabel('Time from cue')
    ylabel('Spike rate (Hz)')
    title('Reward')
    ylim([0 3])
    xlim([-500 2000])
    vline([600], 'b')
    if id == 4
        vline([1100], 'c')
    end
    subplot(2,3,3)
    plot(tt, nanmean(all_earlytrial_rew_lick,2).*(1000./frameRateHz),'k');
    hold on
    plot(tt, nanmean(all_latetrial_rew_lick,2).*(1000./frameRateHz),'b');
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
        shadedErrorBar(tt, nanmean(all_earlytrial_unexp_df,2), (nanstd(all_earlytrial_unexp_df,[],2)./sqrt(totIC)),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_latetrial_unexp_df,2), (nanstd(all_latetrial_unexp_df,[],2)./sqrt(totIC)),'b');
        title('Unexpected')
    else
        shadedErrorBar(tt, nanmean(all_earlytrial_omit_df,2), (nanstd(all_earlytrial_omit_df,[],2)./sqrt(totIC)),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_latetrial_omit_df,2), (nanstd(all_latetrial_omit_df,[],2)./sqrt(totIC)),'b');
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
        shadedErrorBar(tt, nanmean(all_earlytrial_unexp,2).*(1000./frameRateHz), (nanstd(all_earlytrial_unexp,[],2)./sqrt(totIC)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_latetrial_unexp,2).*(1000./frameRateHz), (nanstd(all_latetrial_unexp,[],2)./sqrt(totIC)).*(1000./frameRateHz),'b');
        title('Unexpected')
    else
        shadedErrorBar(tt, nanmean(all_earlytrial_omit,2).*(1000./frameRateHz), (nanstd(all_earlytrial_omit,[],2)./sqrt(totIC)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_latetrial_omit,2).*(1000./frameRateHz), (nanstd(all_latetrial_omit,[],2)./sqrt(totIC)).*(1000./frameRateHz),'b');
        title('Omission')
    end
    xlabel('Time from cue')
    ylabel('Spike rate (Hz)')
    ylim([0 3])
    xlim([-500 2000])
    vline([600], 'b')
    if id == 4
        vline([1100], 'c')
    end
    subplot(2,3,6)
    if id == 3
        plot(tt, nanmean(all_earlytrial_unexp_lick,2).*(1000./frameRateHz), 'k');
        hold on
        plot(tt, nanmean(all_latetrial_unexp_lick,2).*(1000./frameRateHz), 'b');
        title('Unexpected')
    else
        plot(tt, nanmean(all_earlytrial_omit_lick,2).*(1000./frameRateHz), 'k');
        hold on
        plot(tt, nanmean(all_latetrial_omit_lick,2).*(1000./frameRateHz), 'b');
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
    suptitle(['Day ' num2str(id) ': n= ' num2str(nexp) ' mice, ' num2str(totIC) ' dendrites- early (black) vs late (blue) trials']);
    savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_earlyVlateTrial' mouse_str '.fig'])    
    
    if size(all_area_id,2)>10
        totExp = sum(expt_areas,2);
        for i = 1:size(area_list,2)
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
            ylim([0 3])
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
            ylim([0 3])
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
    end
    
    figure;
    subplot(2,2,1)
    shadedErrorBar(tt, nanmean(all_earlylick_rew,2).*(1000./frameRateHz), (nanstd(all_earlylick_rew,[],2)./sqrt(totIC)).*(1000./frameRateHz),'k');
    hold on
    shadedErrorBar(tt, nanmean(all_latelick_rew,2).*(1000./frameRateHz), (nanstd(all_latelick_rew,[],2)./sqrt(totIC)).*(1000./frameRateHz),'b');
    scatter(all_early_rew_time, -0.5.*ones(size(all_early_rew_time)), 'xk');
    scatter(all_late_rew_time, -0.5.*ones(size(all_late_rew_time)), 'xb');
    xlabel('Time from cue')
    ylabel('Spike rate (Hz)')
    title('Rewarded trials')
    ylim([-1 4])
    xlim([-500 2000])
    vline([600], 'b')
    if id == 4
        vline([1100], 'c')
    end
    subplot(2,2,2)
    shadedErrorBar(tt, nanmean(all_earlylick_omit,2).*(1000./frameRateHz), (nanstd(all_earlylick_omit,[],2)./sqrt(totIC)).*(1000./frameRateHz),'k');
    hold on
    shadedErrorBar(tt, nanmean(all_latelick_omit,2).*(1000./frameRateHz), (nanstd(all_latelick_omit,[],2)./sqrt(totIC)).*(1000./frameRateHz),'b');
    scatter(all_early_omit_time, -0.5.*ones(size(all_early_omit_time)), 'xk');
    scatter(all_late_omit_time, -0.5.*ones(size(all_late_omit_time)), 'xb');
    xlabel('Time from cue')
    ylabel('Spike rate (Hz)')
    title('Omission trials')
    ylim([-1 4])
    xlim([-500 2000])
    vline([600], 'b')
    if id == 4
        vline([1100], 'c')
    end
    if doSpr2017sep
        subplot(2,2,3)
        shadedErrorBar(tt, nanmean(all_early_rew_6,2).*(1000./frameRateHz), (nanstd(all_early_rew_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_late_rew_6,2).*(1000./frameRateHz), (nanstd(all_late_rew_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'b');
        scatter(all_early_rew_time_6, -0.5.*ones(size(all_early_rew_time_6)), 'xk');
        scatter(all_late_rew_time_6, -0.5.*ones(size(all_late_rew_time_6)), 'xb');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Spring 2017 mice')
        ylim([-1 4])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end

        subplot(2,2,4)
        shadedErrorBar(tt, nanmean(all_early_omit_6,2).*(1000./frameRateHz), (nanstd(all_early_omit_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_late_omit_6,2).*(1000./frameRateHz), (nanstd(all_late_omit_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'b');
        scatter(all_early_omit_time_6, -0.5.*ones(size(all_early_omit_time_6)), 'xk');
        scatter(all_late_omit_time_6, -0.5.*ones(size(all_late_omit_time_6)), 'xb');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Spring 2017 mice')
        ylim([-1 4])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
    end
    suptitle(['Day ' num2str(id) ': n= ' num2str(nexp) ' mice, ' num2str(totIC) ' dendrites- Early (black) vs late (blue) lick bursts']);
    savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_EarlyVsLateLicks' mouse_str '.fig'])    
    
    if size(all_area_id,2)>10
        totExp = sum(expt_areas,2);
        for i = 1:size(area_list,2)
            ind = find(all_area_id==i);
            totIC_area = length(ind);
            figure;
            subplot(2,2,1)
            shadedErrorBar(tt, nanmean(all_earlylick_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_earlylick_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(all_latelick_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_latelick_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
            scatter(all_early_rew_time(:,find(expt_areas(i,:))), -0.5.*ones(size(all_early_rew_time(:,find(expt_areas(i,:))))), 'xk');
            scatter(all_late_rew_time(:,find(expt_areas(i,:))), -0.5.*ones(size(all_late_rew_time(:,find(expt_areas(i,:))))), 'xb');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title('Rewarded trials')
            ylim([-1 4])
            xlim([-500 2000])
            vline([600], 'b')
            if id == 4
                vline([1100], 'c')
            end
            subplot(2,2,2)
            shadedErrorBar(tt, nanmean(all_earlylick_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_earlylick_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(all_latelick_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_latelick_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
            scatter(all_early_omit_time(:,find(expt_areas(i,:))), -0.5.*ones(size(all_early_omit_time(:,find(expt_areas(i,:))))), 'xk');
            scatter(all_late_omit_time(:,find(expt_areas(i,:))), -0.5.*ones(size(all_late_omit_time(:,find(expt_areas(i,:))))), 'xb');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title('Omission trials')
            ylim([-1 4])
            xlim([-500 2000])
            vline([600], 'b')
            if id == 4
                vline([1100], 'c')
            end
            suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice, ' num2str(totIC_area) ' dendrites- Early (black) vs late (blue) lick bursts']);
            savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_EarlyVsLateLicks_' mouse_str '.fig'])     
        end
    end
    
    figure;
    if size(all_short_omit,2)
        subplot(3,1,1)
        shadedErrorBar(tt, nanmean(all_short_omit,2).*(1000./frameRateHz), (nanstd(all_short_omit,[],2)./sqrt(totIC)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_long_omit,2).*(1000./frameRateHz), (nanstd(all_long_omit,[],2)./sqrt(totIC)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([-1 4])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        title('Omits after short (black) and long (blue) interval')
    end
    if size(all_short_unexp,2)
        subplot(3,1,2)
        shadedErrorBar(tt, nanmean(all_short_unexp,2).*(1000./frameRateHz), (nanstd(all_short_unexp,[],2)./sqrt(totIC)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_long_unexp,2).*(1000./frameRateHz), (nanstd(all_long_unexp,[],2)./sqrt(totIC)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([-1 4])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        title('Unexpected reward after short (black) and long (blue) interval')
    end
    if size(all_preomit,2)
        subplot(3,1,3)
        shadedErrorBar(tt, nanmean(all_preomit,2).*(1000./frameRateHz), (nanstd(all_preomit,[],2)./sqrt(totIC)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_postomit,2).*(1000./frameRateHz), (nanstd(all_postomit,[],2)./sqrt(totIC)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        ylim([-1 4])
        xlim([-500 2000])
        vline([600], 'b')
        if id == 4
            vline([1100], 'c')
        end
        title('Reward before (black) and after (blue) omit trial')
    end
    suptitle(['Day ' num2str(id) ': n= ' num2str(nexp) ' mice, ' num2str(totIC) ' dendrites']);
    savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_trialByTrialAnalysis' mouse_str '.fig'])    
    
    if doSpr2017sep
        figure;
        if size(all_short_omit_6,2)
            subplot(3,1,1)
            shadedErrorBar(tt, nanmean(all_short_omit_6,2).*(1000./frameRateHz), (nanstd(all_short_omit_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(all_long_omit_6,2).*(1000./frameRateHz), (nanstd(all_long_omit_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([-1 4])
            xlim([-500 2000])
            vline([600], 'b')
            if id == 4
                vline([1100], 'c')
            end
            title('Omits after short (black) and long (blue) interval')
        end
        if size(all_short_unexp_6,2)
            subplot(3,1,2)
            shadedErrorBar(tt, nanmean(all_short_unexp_6,2).*(1000./frameRateHz), (nanstd(all_short_unexp_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(all_long_unexp_6,2).*(1000./frameRateHz), (nanstd(all_long_unexp_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([-1 4])
            xlim([-500 2000])
            vline([600], 'b')
            if id == 4
                vline([1100], 'c')
            end
            title('Unexpected reward after short (black) and long (blue) interval')
        end
        if size(all_preomit_6,2)
            subplot(3,1,3)
            shadedErrorBar(tt, nanmean(all_preomit_6,2).*(1000./frameRateHz), (nanstd(all_preomit_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'k');
            hold on
            shadedErrorBar(tt, nanmean(all_postomit_6,2).*(1000./frameRateHz), (nanstd(all_postomit_6,[],2)./sqrt(totIC_6)).*(1000./frameRateHz),'b');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            ylim([-1 4])
            xlim([-500 2000])
            vline([600], 'b')
            if id == 4
                vline([1100], 'c')
            end
            title('Reward before (black) and after (blue) omit trial')
        end
        suptitle(['Day ' num2str(id) ': n= 6 mice, ' num2str(totIC_6) ' dendrites- Spring 2017 mice']);
        savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_trialByTrialAnalysis_spring2017mice' mouse_str_6 '.fig'])    
    end
    
    figure;
    hist(all_pct_precue_burst,[0:.01:.5])
    xlabel('Fract trials with pre-reward lick bursts')
    ylabel('Mice')
    title(['Day ' num2str(id) ': n= ' num2str(nexp) ' mice- Pre-reward lick burst frequency']);
    savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_preRewardLicks' mouse_str '.fig'])    
    
    prerew_resp_cell = [];
    postrew_resp_cell = [];
    prerew_frames = round(400./frameRateHz);
    rewdelay_frames = round(600./frameRateHz);
    base_avg = mean(all_rew(1:prewin_frames,:),1);
    base_std = std(all_rew(1:prewin_frames,:),[],1);
    for iC = 1:totIC
        prediff_ind = diff(intersect(prewin_frames+1:prewin_frames+prerew_frames, find(all_rew(:,iC)>base_avg+(base_std.*1))))';
        if length(prediff_ind)
            if strfind(prediff_ind, [1 1])
                prerew_resp_cell = [prerew_resp_cell iC];
            end
        end
        postdiff_ind = diff(intersect(prewin_frames+rewdelay_frames+1:prewin_frames+prerew_frames+rewdelay_frames, find(all_rew(:,iC)>base_avg+(base_std.*1))))';
        if length(postdiff_ind)
            if find(postdiff_ind==1)
                postrew_resp_cell = [postrew_resp_cell iC];
            end
        end
    end
    prerew_notresp_cell = setdiff(1:totIC,prerew_resp_cell);
    postrew_notresp_cell = setdiff(1:totIC,postrew_resp_cell);
    
    figure; 
    subplot(2,2,1)
    shadedErrorBar(tt, nanmean(all_rew(:,prerew_resp_cell),2).*(1000./frameRateHz), (nanstd(all_rew(:,prerew_resp_cell),[],2)./sqrt(length(prerew_resp_cell))).*(1000./frameRateHz),'k');
    hold on
    if id == 3
        shadedErrorBar(tt, nanmean(all_unexp(:,prerew_resp_cell),2).*(1000./frameRateHz), (nanstd(all_unexp(:,prerew_resp_cell),[],2)./sqrt(length(prerew_resp_cell))).*(1000./frameRateHz),'g');
    else
        shadedErrorBar(tt, nanmean(all_omit(:,prerew_resp_cell),2).*(1000./frameRateHz), (nanstd(all_omit(:,prerew_resp_cell),[],2)./sqrt(length(prerew_resp_cell))).*(1000./frameRateHz),'r');
    end
    xlabel('Time from cue')
    ylabel('Spike rate (Hz)')
    title(['Pre-reward responsive: n = ' num2str(length(prerew_resp_cell))])
    ylim([-1 4])
    vline([600], 'b')
    if id == 4
        vline([1100], 'c')
    end
    xlim([-500 2000])
    subplot(2,2,3)
    shadedErrorBar(tt, nanmean(all_rew(:,prerew_notresp_cell),2).*(1000./frameRateHz), (nanstd(all_rew(:,prerew_notresp_cell),[],2)./sqrt(length(prerew_notresp_cell))).*(1000./frameRateHz),'k');
    hold on
    if id == 3
        shadedErrorBar(tt, nanmean(all_unexp(:,prerew_notresp_cell),2).*(1000./frameRateHz), (nanstd(all_unexp(:,prerew_notresp_cell),[],2)./sqrt(length(prerew_notresp_cell))).*(1000./frameRateHz),'g');
    else
        shadedErrorBar(tt, nanmean(all_omit(:,prerew_notresp_cell),2).*(1000./frameRateHz), (nanstd(all_omit(:,prerew_notresp_cell),[],2)./sqrt(length(prerew_notresp_cell))).*(1000./frameRateHz),'r');
    end
    xlabel('Time from cue')
    ylabel('Spike rate (Hz)')
    title(['Not pre-reward responsive: n = ' num2str(length(prerew_notresp_cell))])
    ylim([-1 4])
    xlim([-500 2000])
    vline([600], 'b')
    if id == 4
        vline([1100], 'c')
    end
    subplot(2,2,2)
    shadedErrorBar(tt, nanmean(all_rew(:,postrew_resp_cell),2).*(1000./frameRateHz), (nanstd(all_rew(:,postrew_resp_cell),[],2)./sqrt(length(postrew_resp_cell))).*(1000./frameRateHz),'k');
    hold on
    if id == 3
        shadedErrorBar(tt, nanmean(all_unexp(:,postrew_resp_cell),2).*(1000./frameRateHz), (nanstd(all_unexp(:,postrew_resp_cell),[],2)./sqrt(length(postrew_resp_cell))).*(1000./frameRateHz),'g');
    else
        shadedErrorBar(tt, nanmean(all_omit(:,postrew_resp_cell),2).*(1000./frameRateHz), (nanstd(all_omit(:,postrew_resp_cell),[],2)./sqrt(length(postrew_resp_cell))).*(1000./frameRateHz),'r');
    end
    xlabel('Time from cue')
    ylabel('Spike rate (Hz)')
    title(['Post-reward responsive: n = ' num2str(length(postrew_resp_cell))])
    ylim([-1 4])
    vline([600], 'b')
    if id == 4
        vline([1100], 'c')
    end
    xlim([-500 2000])
    subplot(2,2,4)
    shadedErrorBar(tt, nanmean(all_rew(:,postrew_notresp_cell),2).*(1000./frameRateHz), (nanstd(all_rew(:,postrew_notresp_cell),[],2)./sqrt(length(postrew_notresp_cell))).*(1000./frameRateHz),'k');
    hold on
    if id == 3
        shadedErrorBar(tt, nanmean(all_unexp(:,postrew_notresp_cell),2).*(1000./frameRateHz), (nanstd(all_unexp(:,postrew_notresp_cell),[],2)./sqrt(length(postrew_notresp_cell))).*(1000./frameRateHz),'g');
    else
        shadedErrorBar(tt, nanmean(all_omit(:,postrew_notresp_cell),2).*(1000./frameRateHz), (nanstd(all_omit(:,postrew_notresp_cell),[],2)./sqrt(length(postrew_notresp_cell))).*(1000./frameRateHz),'r');
    end
    xlabel('Time from cue')
    ylabel('Spike rate (Hz)')
    title(['Not post-reward responsive: n = ' num2str(length(postrew_notresp_cell))])
    ylim([-1 4])
    xlim([-500 2000])
    vline([600], 'b')
    if id == 4
        vline([1100], 'c')
    end
    if id == 3
        suptitle(['Day ' num2str(id) ': n= ' num2str(nexp) ' mice- Reward (black); unexpected (green)'])
    else
        suptitle(['Day ' num2str(id) ': n= ' num2str(nexp) ' mice- Reward (black); omit (red)'])
    end
    savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_preVpostRespCells' mouse_str '.fig'])
    
    figure; 
    subplot(2,3,1)
    shadedErrorBar(tl, nanmean(all_postrew_lick_rew,2).*(1000./frameRateHz), (nanstd(all_postrew_lick_rew,[],2).*(1000./frameRateHz))./sqrt(size(all_postrew_lick_rew,2)),'k');
    hold on
    if length(all_postrew_lick_omit)>10
        shadedErrorBar(tl, nanmean(all_postrew_lick_omit,2).*(1000./frameRateHz), (nanstd(all_postrew_lick_omit,[],2).*(1000./frameRateHz))./sqrt(size(all_postrew_lick_omit,2)),'r');
    end
    xlabel('Time from lick')
    ylabel('Spike rate (Hz)')
    ylim([0 4])
    title(['Reward (black, n = ' num2str(length(all_postrew_lick_rew)) '); omit (red, n = ' num2str(length(all_postrew_lick_omit)) ')'])
    subplot(2,3,2)
    shadedErrorBar(tl, nanmean(all_early_postrew_lick_rew,2).*(1000./frameRateHz), (nanstd(all_early_postrew_lick_rew,[],2).*(1000./frameRateHz))./sqrt(size(all_early_postrew_lick_rew,2)),'k');
    hold on
    shadedErrorBar(tl, nanmean(all_late_postrew_lick_rew,2).*(1000./frameRateHz), (nanstd(all_late_postrew_lick_rew,[],2).*(1000./frameRateHz))./sqrt(size(all_late_postrew_lick_rew,2)),'b');
    xlabel('Time from lick')
    ylabel('Spike rate (Hz)')
    ylim([0 4])
    title(['Early licks (black); late licks (blue)'])
    subplot(2,3,4)
    shadedErrorBar(tl, nanmean(all_postrew_lick_rew_lick,2).*(1000./frameRateHz), (nanstd(all_postrew_lick_rew_lick,[],2).*(1000./frameRateHz))./sqrt(size(all_postrew_lick_rew_lick,2)),'k');
    hold on
    if length(all_postrew_lick_omit)>10
        shadedErrorBar(tl, nanmean(all_postrew_lick_omit_lick,2).*(1000./frameRateHz), (nanstd(all_postrew_lick_omit_lick,[],2).*(1000./frameRateHz))./sqrt(size(all_postrew_lick_omit_lick,2)),'r');
    end
    xlabel('Time from lick')
    ylabel('Lick rate (Hz)')
    ylim([0 inf])
    subplot(2,3,5)
    shadedErrorBar(tl, nanmean(all_early_postrew_lick_rew_lick,2).*(1000./frameRateHz), (nanstd(all_early_postrew_lick_rew_lick,[],2).*(1000./frameRateHz))./sqrt(size(all_early_postrew_lick_rew_lick,2)),'k');
    hold on
    shadedErrorBar(tl, nanmean(all_late_postrew_lick_rew_lick,2).*(1000./frameRateHz), (nanstd(all_late_postrew_lick_rew_lick,[],2).*(1000./frameRateHz))./sqrt(size(all_late_postrew_lick_rew_lick,2)),'b');
    xlabel('Time from lick')
    ylabel('Lick rate (Hz)')
    ylim([0 inf])
    if id == 3
        subplot(2,3,3)
        shadedErrorBar(tl, nanmean(all_postrew_lick_unexp,2).*(1000./frameRateHz), (nanstd(all_postrew_lick_unexp,[],2).*(1000./frameRateHz))./sqrt(size(all_postrew_lick_unexp,2)),'g');
        xlabel('Time from lick')
        ylabel('Spike rate (Hz)')
        ylim([0 4])
        title(['Unexpected reward'])
        subplot(2,3,6)
        shadedErrorBar(tl, nanmean(all_postrew_lick_unexp_lick,2).*(1000./frameRateHz), (nanstd(all_postrew_lick_unexp_lick,[],2).*(1000./frameRateHz))./sqrt(size(all_postrew_lick_unexp_lick,2)),'g');
        xlabel('Time from lick')
        ylabel('Lick rate (Hz)')
        ylim([0 inf])
    end
    suptitle(['Post-reward lick burst aligned spiking'])
    savefig(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_postRewLickAlignSpiking' mouse_str '.fig'])

end