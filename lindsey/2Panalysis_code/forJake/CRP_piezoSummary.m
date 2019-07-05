clear all
close all
CRP_expt_list_LS2019
%area_list = {'C1','C2'};
area_list = {'LS'};
share_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\public\ClassicConditioningPaper\LS2019';
if ~exist(share_out)
    mkdir(share_out)
end
for id = [1:2 5]
    fprintf(['Day: ' num2str(id) '\n'])
    lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';
    jake_dir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\jake\Analysis\Cue_reward_pairing_analysis\CC_analysis_2P';
    nexp = size(expt(id).date,1);
    mouse_str = [];
    all_area_id = [];
    expt_areas = zeros(length(area_list),nexp);
    all_rew = [];
    all_omit = [];
    all_rew_piezo = [];
    all_omit_piezo = [];
    all_low25piezo_prerew = [];
    all_high25piezo_prerew = [];
    all_low25piezo_postrew = [];
    all_high25piezo_postrew = [];
    all_low25piezo_preomit = [];
    all_high25piezo_preomit = [];
    all_low25piezo_postomit = [];
    all_high25piezo_postomit = [];
    all_low10piezo_prerew = [];
    all_high10piezo_prerew = [];
    all_low10piezo_postrew = [];
    all_high10piezo_postrew = [];
    all_low10piezo_preomit = [];
    all_high10piezo_preomit = [];
    all_low10piezo_postomit = [];
    all_high10piezo_postomit = [];
    all_earlypiezo_rew = [];
    all_latepiezo_rew = [];
    all_earlypiezo_omit = [];
    all_latepiezo_omit = [];
    all_earlypiezo_rew_piezo = [];
    all_latepiezo_rew_piezo = [];
    all_earlypiezo_omit_piezo = [];
    all_latepiezo_omit_piezo = [];
    all_allearlypiezo25_rew = [];
    all_alllatepiezo25_rew = [];
    all_allearlypiezo25_omit = [];
    all_alllatepiezo25_omit = [];
    all_allearlypiezo25_rew_piezo = [];
    all_alllatepiezo25_rew_piezo = [];
    all_allearlypiezo25_omit_piezo = [];
    all_alllatepiezo25_omit_piezo = [];
    all_allearlypiezo10_rew = [];
    all_alllatepiezo10_rew = [];
    all_allearlypiezo10_omit = [];
    all_alllatepiezo10_omit = [];
    all_allearlypiezo10_rew_piezo = [];
    all_alllatepiezo10_rew_piezo = [];
    all_allearlypiezo10_omit_piezo = [];
    all_alllatepiezo10_omit_piezo = [];
    all_high25piezo_prerew_piezo = [];
    all_low25piezo_prerew_piezo = [];
    all_high25piezo_postrew_piezo = [];
    all_low25piezo_postrew_piezo = [];
    all_high25piezo_preomit_piezo = [];
    all_low25piezo_preomit_piezo = [];
    all_high25piezo_postomit_piezo = [];
    all_low25piezo_postomit_piezo = [];
    all_high10piezo_prerew_piezo = [];
    all_low10piezo_prerew_piezo = [];
    all_high10piezo_postrew_piezo = [];
    all_low10piezo_postrew_piezo = [];
    all_high10piezo_preomit_piezo = [];
    all_low10piezo_preomit_piezo = [];
    all_high10piezo_postomit_piezo = [];
    all_low10piezo_postomit_piezo = [];

    for iexp = 1:nexp
        mouse = expt(id).mouse(iexp,:);
        if find(mouse == ' ')
            ind = find(mouse == ' ');
            mouse(ind) = [];
        end
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];
        mouse_str = [mouse_str '_' mouse];
        
        load(fullfile(lg_out,img_fn, [img_fn '_targetAlign.mat']))
        load(fullfile(lg_out,img_fn, [img_fn '_cueAlignPiezo.mat']))
        if exist(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']))
            load(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']))
            area_temp = zeros(1,size(targetAligndFoverF,2));
            for i = 1:2
                if ~isnan(cell2mat(expt(id).areas{iexp}(i)))
                    x = find(strcmp(area_list, expt(id).areas{iexp}(i)));
                    ind = find(maskCat==i);
                    area_temp(1,ind) = x;
                    expt_areas(x,iexp) = 1;
                else
                    ind = find(maskCat==i);
                    area_temp(1,ind) = NaN;
                end    
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
        
        rew_piezo_avg =nanmean(targetAlign_piezo(:,ind_rew),2);
        omit_piezo_avg =nanmean(targetAlign_piezo(:,ind_omit),2);

        earlypiezo_rew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_earlypiezo_rew)),3);
        latepiezo_rew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_latepiezo_rew)),3);
        earlypiezo_omit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_earlypiezo_omit)),3);
        latepiezo_omit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_latepiezo_omit)),3);
        
        earlypiezo_rew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_earlypiezo_rew)),2);
        latepiezo_rew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_latepiezo_rew)),2);
        earlypiezo_omit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_earlypiezo_omit)),2);
        latepiezo_omit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_latepiezo_omit)),2);
        
        allearlypiezo25_rew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_allearlypiezo25_rew)),3);
        alllatepiezo25_rew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_alllatepiezo25_rew)),3);
        allearlypiezo25_omit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_allearlypiezo25_omit)),3);
        alllatepiezo25_omit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_alllatepiezo25_omit)),3);
        
        allearlypiezo25_rew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_allearlypiezo25_rew)),2);
        alllatepiezo25_rew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_alllatepiezo25_rew)),2);
        allearlypiezo25_omit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_allearlypiezo25_omit)),2);
        alllatepiezo25_omit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_alllatepiezo25_omit)),2);
        
        allearlypiezo10_rew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_allearlypiezo10_rew)),3);
        alllatepiezo10_rew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_alllatepiezo10_rew)),3);
        allearlypiezo10_omit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_allearlypiezo10_omit)),3);
        alllatepiezo10_omit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_alllatepiezo10_omit)),3);
        
        allearlypiezo10_rew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_allearlypiezo10_rew)),2);
        alllatepiezo10_rew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_alllatepiezo10_rew)),2);
        allearlypiezo10_omit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_allearlypiezo10_omit)),2);
        alllatepiezo10_omit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_alllatepiezo10_omit)),2);
        
        low25piezo_prerew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_low25piezo_prerew)),3);
        high25piezo_prerew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_high25piezo_prerew)),3);
        low25piezo_postrew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_low25piezo_postrew)),3);
        high25piezo_postrew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_high25piezo_postrew)),3);
        low25piezo_preomit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_low25piezo_preomit)),3);
        high25piezo_preomit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_high25piezo_preomit)),3);
        low25piezo_postomit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_low25piezo_postomit)),3);
        high25piezo_postomit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_high25piezo_postomit)),3);
        
        low25piezo_prerew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_low25piezo_prerew)),2);
        high25piezo_prerew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_high25piezo_prerew)),2);
        low25piezo_postrew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_low25piezo_postrew)),2);
        high25piezo_postrew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_high25piezo_postrew)),2);
        low25piezo_preomit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_low25piezo_preomit)),2);
        high25piezo_preomit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_high25piezo_preomit)),2);
        low25piezo_postomit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_low25piezo_postomit)),2);
        high25piezo_postomit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_high25piezo_postomit)),2);
        
        low10piezo_prerew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_low10piezo_prerew)),3);
        high10piezo_prerew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_high10piezo_prerew)),3);
        low10piezo_postrew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_low10piezo_postrew)),3);
        high10piezo_postrew_avg = nanmean(targetAlign_events(:,:,ind_rew(ind_high10piezo_postrew)),3);
        low10piezo_preomit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_low10piezo_preomit)),3);
        high10piezo_preomit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_high10piezo_preomit)),3);
        low10piezo_postomit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_low10piezo_postomit)),3);
        high10piezo_postomit_avg = nanmean(targetAlign_events(:,:,ind_omit(ind_high10piezo_postomit)),3);
        
        low10piezo_prerew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_low10piezo_prerew)),2);
        high10piezo_prerew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_high10piezo_prerew)),2);
        low10piezo_postrew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_low10piezo_postrew)),2);
        high10piezo_postrew_piezoavg = nanmean(targetAlign_piezo(:,ind_rew(ind_high10piezo_postrew)),2);
        low10piezo_preomit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_low10piezo_preomit)),2);
        high10piezo_preomit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_high10piezo_preomit)),2);
        low10piezo_postomit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_low10piezo_postomit)),2);
        high10piezo_postomit_piezoavg = nanmean(targetAlign_piezo(:,ind_omit(ind_high10piezo_postomit)),2);
        
        n_rew = size(ind_rew,2);
        n_omit = size(ind_omit,2);
        n_unexp = size(ind_unexp,2);
        nIC = size(rew_avg,2);
                
        if ~exist([share_out '\CC_summary\Day' num2str(id)])
            mkdir([share_out '\CC_summary\Day' num2str(id)])
        end

        all_rew = [all_rew rew_avg];
        all_omit = [all_omit omit_avg];
        
        all_rew_piezo = [all_rew_piezo rew_piezo_avg];
        all_omit_piezo = [all_omit_piezo omit_piezo_avg];
        
        all_earlypiezo_rew = [all_earlypiezo_rew earlypiezo_rew_avg];
        all_latepiezo_rew = [all_latepiezo_rew latepiezo_rew_avg];
        all_earlypiezo_omit = [all_earlypiezo_omit earlypiezo_omit_avg];
        all_latepiezo_omit = [all_latepiezo_omit latepiezo_omit_avg];
        all_allearlypiezo25_rew = [all_allearlypiezo25_rew allearlypiezo25_rew_avg];
        all_alllatepiezo25_rew = [all_alllatepiezo25_rew alllatepiezo25_rew_avg];
        all_allearlypiezo25_omit = [all_allearlypiezo25_omit allearlypiezo25_omit_avg];
        all_alllatepiezo25_omit = [all_alllatepiezo25_omit alllatepiezo25_omit_avg];
        all_allearlypiezo10_rew = [all_allearlypiezo10_rew allearlypiezo10_rew_avg];
        all_alllatepiezo10_rew = [all_alllatepiezo10_rew alllatepiezo10_rew_avg];
        all_allearlypiezo10_omit = [all_allearlypiezo10_omit allearlypiezo10_omit_avg];
        all_alllatepiezo10_omit = [all_alllatepiezo10_omit alllatepiezo10_omit_avg];
        
        all_low25piezo_prerew = [all_low25piezo_prerew low25piezo_prerew_avg];
        all_high25piezo_prerew = [all_high25piezo_prerew high25piezo_prerew_avg];
        all_low25piezo_postrew = [all_low25piezo_postrew low25piezo_postrew_avg];
        all_high25piezo_postrew = [all_high25piezo_postrew high25piezo_postrew_avg];
        all_low25piezo_preomit = [all_low25piezo_preomit low25piezo_preomit_avg];
        all_high25piezo_preomit = [all_high25piezo_preomit high25piezo_preomit_avg];
        all_low25piezo_postomit = [all_low25piezo_postomit low25piezo_postomit_avg];
        all_high25piezo_postomit = [all_high25piezo_postomit high25piezo_postomit_avg];
        
        all_earlypiezo_rew_piezo = [all_earlypiezo_rew_piezo earlypiezo_rew_piezoavg];
        all_latepiezo_rew_piezo = [all_latepiezo_rew_piezo latepiezo_rew_piezoavg];
        all_earlypiezo_omit_piezo = [all_earlypiezo_omit_piezo earlypiezo_omit_piezoavg];
        all_latepiezo_omit_piezo = [all_latepiezo_omit_piezo latepiezo_omit_piezoavg];
        
        all_high25piezo_prerew_piezo = [all_high25piezo_prerew_piezo high25piezo_prerew_piezoavg];
        all_low25piezo_prerew_piezo = [all_low25piezo_prerew_piezo low25piezo_prerew_piezoavg];
        all_high25piezo_postrew_piezo = [all_high25piezo_postrew_piezo high25piezo_postrew_piezoavg];
        all_low25piezo_postrew_piezo = [all_low25piezo_postrew_piezo low25piezo_postrew_piezoavg];
        all_high25piezo_preomit_piezo = [all_high25piezo_preomit_piezo high25piezo_preomit_piezoavg];
        all_low25piezo_preomit_piezo = [all_low25piezo_preomit_piezo low25piezo_preomit_piezoavg];
        all_high25piezo_postomit_piezo = [all_high25piezo_postomit_piezo high25piezo_postomit_piezoavg];
        all_low25piezo_postomit_piezo = [all_low25piezo_postomit_piezo low25piezo_postomit_piezoavg];
        
        all_low10piezo_prerew = [all_low10piezo_prerew low10piezo_prerew_avg];
        all_high10piezo_prerew = [all_high10piezo_prerew high10piezo_prerew_avg];
        all_low10piezo_postrew = [all_low10piezo_postrew low10piezo_postrew_avg];
        all_high10piezo_postrew = [all_high10piezo_postrew high10piezo_postrew_avg];
        all_low10piezo_preomit = [all_low10piezo_preomit low10piezo_preomit_avg];
        all_high10piezo_preomit = [all_high10piezo_preomit high10piezo_preomit_avg];
        all_low10piezo_postomit = [all_low10piezo_postomit low10piezo_postomit_avg];
        all_high10piezo_postomit = [all_high10piezo_postomit high10piezo_postomit_avg];
        
        all_allearlypiezo25_rew_piezo = [all_allearlypiezo25_rew_piezo allearlypiezo25_rew_piezoavg];
        all_alllatepiezo25_rew_piezo = [all_alllatepiezo25_rew_piezo alllatepiezo25_rew_piezoavg];
        all_allearlypiezo25_omit_piezo = [all_allearlypiezo25_omit_piezo allearlypiezo25_omit_piezoavg];
        all_alllatepiezo25_omit_piezo = [all_alllatepiezo25_omit_piezo alllatepiezo25_omit_piezoavg];
        all_allearlypiezo10_rew_piezo = [all_allearlypiezo10_rew_piezo allearlypiezo10_rew_piezoavg];
        all_alllatepiezo10_rew_piezo = [all_alllatepiezo10_rew_piezo alllatepiezo10_rew_piezoavg];
        all_allearlypiezo10_omit_piezo = [all_allearlypiezo10_omit_piezo allearlypiezo10_omit_piezoavg];
        all_alllatepiezo10_omit_piezo = [all_alllatepiezo10_omit_piezo alllatepiezo10_omit_piezoavg];
        
        all_high10piezo_prerew_piezo = [all_high10piezo_prerew_piezo high10piezo_prerew_piezoavg];
        all_low10piezo_prerew_piezo = [all_low10piezo_prerew_piezo low10piezo_prerew_piezoavg];
        all_high10piezo_postrew_piezo = [all_high10piezo_postrew_piezo high10piezo_postrew_piezoavg];
        all_low10piezo_postrew_piezo = [all_low10piezo_postrew_piezo low10piezo_postrew_piezoavg];
        all_high10piezo_preomit_piezo = [all_high10piezo_preomit_piezo high10piezo_preomit_piezoavg];
        all_low10piezo_preomit_piezo = [all_low10piezo_preomit_piezo low10piezo_preomit_piezoavg];
        all_high10piezo_postomit_piezo = [all_high10piezo_postomit_piezo high10piezo_postomit_piezoavg];
        all_low10piezo_postomit_piezo = [all_low10piezo_postomit_piezo low10piezo_postomit_piezoavg];
        
        all_HL_piezo(iexp) = HL_piezo;
    end
        
    totIC = size(all_area_id,2);
    totExp = sum(expt_areas,2);
    all_HL_piezo = concatenateStructuresLG(all_HL_piezo);

    for i = find(sum(expt_areas,2)')
        ind = find(all_area_id==i);
        totIC_area = length(ind);
        figure;
        shadedErrorBar(tt, nanmean(all_rew_piezo(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), (nanstd(all_rew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_omit_piezo(:,find(expt_areas(i,:))),2).*(1000./frameRateHz), (nanstd(all_omit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))).*(1000./frameRateHz),'r');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice- Reward (black) vs omit (red)']);
        savefig([share_out '\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_RewVsOmitPiezo_' mouse_str '.fig'])     
        
        figure;
        subplot(2,2,1)
        shadedErrorBar(tt, nanmean(all_earlypiezo_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_earlypiezo_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_latepiezo_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_latepiezo_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        hold on
        errorbarxy(nanmean(all_HL_piezo.early_rew(:,find(expt_areas(i,:))).*(1000./frameRateHz),2), 0, nanstd(all_HL_piezo.early_rew(:,find(expt_areas(i,:))).*(1000./frameRateHz),[],2)./sqrt(sum(expt_areas(i,:),2)),0, {'ok', 'k', 'k'});
        hold on
        errorbarxy(nanmean(all_HL_piezo.late_rew(:,find(expt_areas(i,:))).*(1000./frameRateHz),2), 0, nanstd(all_HL_piezo.late_rew(:,find(expt_areas(i,:))).*(1000./frameRateHz),[],2)./sqrt(sum(expt_areas(i,:),2)),0, {'ob', 'b', 'b'});
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Rewarded trials')
        ylim([-1 12])
        vline(600)
        
        subplot(2,2,3)
        shadedErrorBar(tt, nanmean(all_earlypiezo_rew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_earlypiezo_rew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_latepiezo_rew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_latepiezo_rew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        
        subplot(2,2,2)
        shadedErrorBar(tt, nanmean(all_earlypiezo_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_earlypiezo_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_latepiezo_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_latepiezo_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        hold on
        errorbarxy(nanmean(all_HL_piezo.early_omit(:,find(expt_areas(i,:))).*(1000./frameRateHz),2),0, nanstd(all_HL_piezo.early_omit(:,find(expt_areas(i,:))).*(1000./frameRateHz),[],2)./sqrt(sum(expt_areas(i,:),2)),0, {'ok', 'k', 'k'});
        hold on
        errorbarxy(nanmean(all_HL_piezo.late_omit(:,find(expt_areas(i,:))).*(1000./frameRateHz),2), 0, nanstd(all_HL_piezo.late_omit(:,find(expt_areas(i,:))).*(1000./frameRateHz),[],2)./sqrt(sum(expt_areas(i,:),2)),0, {'ob', 'b', 'b'});
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omission trials')
        ylim([-1 12])
        
        subplot(2,2,4)
        shadedErrorBar(tt, nanmean(all_earlypiezo_omit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_earlypiezo_omit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_latepiezo_omit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_latepiezo_omit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice, ' num2str(totIC_area) ' dendrites- Early (black) vs late (blue) movements']);
        savefig([share_out '\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_EarlyVsLatePiezo_' mouse_str '.fig'])     
    
        figure;
        subplot(2,2,1)
        shadedErrorBar(tt, nanmean(all_allearlypiezo25_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_allearlypiezo25_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_alllatepiezo25_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_alllatepiezo25_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        ylabel('Spike rate (Hz)')
        title('Rewarded trials')
        ylim([-1 12])
        vline(600)
        
        subplot(2,2,3)
        shadedErrorBar(tt, nanmean(all_allearlypiezo25_rew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_allearlypiezo25_rew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_alllatepiezo25_rew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_alllatepiezo25_rew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        
        subplot(2,2,2)
        shadedErrorBar(tt, nanmean(all_allearlypiezo25_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_allearlypiezo25_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_alllatepiezo25_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_alllatepiezo25_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omission trials')
        ylim([-1 12])
        
        subplot(2,2,4)
        shadedErrorBar(tt, nanmean(all_allearlypiezo25_omit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_allearlypiezo25_omit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_alllatepiezo25_omit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_alllatepiezo25_omit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice, ' num2str(totIC_area) ' dendrites- All Early25 (black) vs late25 (blue) movements']);
        savefig([share_out '\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_AllEarlyVsLatePiezo25_' mouse_str '.fig'])     
    
        figure;
        subplot(2,2,1)
        shadedErrorBar(tt, nanmean(all_allearlypiezo10_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_allearlypiezo10_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_alllatepiezo10_rew(:,ind),2).*(1000./frameRateHz), (nanstd(all_alllatepiezo10_rew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        ylabel('Spike rate (Hz)')
        title('Rewarded trials')
        ylim([-1 12])
        vline(600)
        
        subplot(2,2,3)
        shadedErrorBar(tt, nanmean(all_allearlypiezo10_rew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_allearlypiezo10_rew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_alllatepiezo10_rew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_alllatepiezo10_rew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        
        subplot(2,2,2)
        shadedErrorBar(tt, nanmean(all_allearlypiezo10_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_allearlypiezo10_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_alllatepiezo10_omit(:,ind),2).*(1000./frameRateHz), (nanstd(all_alllatepiezo10_omit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omission trials')
        ylim([-1 12])
        
        subplot(2,2,4)
        shadedErrorBar(tt, nanmean(all_allearlypiezo10_omit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_allearlypiezo10_omit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_alllatepiezo10_omit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_alllatepiezo10_omit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice, ' num2str(totIC_area) ' dendrites- All Early10 (black) vs late10 (blue) movements']);
        savefig([share_out '\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_AllEarlyVsLatePiezo10_' mouse_str '.fig'])     
    
        figure;
        subplot(2,2,1)
        shadedErrorBar(tt, nanmean(all_high25piezo_prerew(:,ind),2).*(1000./frameRateHz), (nanstd(all_high25piezo_prerew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low25piezo_prerew(:,ind),2).*(1000./frameRateHz), (nanstd(all_low25piezo_prerew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Rewarded trials')
        ylim([-1 12])
        vline(600)
        
        subplot(2,2,3)
        shadedErrorBar(tt, nanmean(all_high25piezo_prerew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_high25piezo_prerew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low25piezo_prerew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_low25piezo_prerew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        
        subplot(2,2,2)
        shadedErrorBar(tt, nanmean(all_high25piezo_preomit(:,ind),2).*(1000./frameRateHz), (nanstd(all_high25piezo_preomit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low25piezo_preomit(:,ind),2).*(1000./frameRateHz), (nanstd(all_low25piezo_preomit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omission trials')
        ylim([-1 12])
        
        subplot(2,2,4)
        shadedErrorBar(tt, nanmean(all_high25piezo_preomit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_high25piezo_preomit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low25piezo_preomit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_low25piezo_preomit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice, ' num2str(totIC_area) ' dendrites- High25 (black) vs low25 (blue) pre-reward movements']);
        savefig([share_out '\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_HighVsLow25PreRewPiezo_' mouse_str '.fig'])     
        
        figure;
        subplot(2,2,1)
        shadedErrorBar(tt, nanmean(all_high10piezo_prerew(:,ind),2).*(1000./frameRateHz), (nanstd(all_high10piezo_prerew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low10piezo_prerew(:,ind),2).*(1000./frameRateHz), (nanstd(all_low10piezo_prerew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Rewarded trials')
        ylim([-1 12])
        vline(600)
        
        subplot(2,2,3)
        shadedErrorBar(tt, nanmean(all_high10piezo_prerew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_high10piezo_prerew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low10piezo_prerew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_low10piezo_prerew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        
        subplot(2,2,2)
        shadedErrorBar(tt, nanmean(all_high10piezo_preomit(:,ind),2).*(1000./frameRateHz), (nanstd(all_high10piezo_preomit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low10piezo_preomit(:,ind),2).*(1000./frameRateHz), (nanstd(all_low10piezo_preomit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omission trials')
        ylim([-1 12])
        
        subplot(2,2,4)
        shadedErrorBar(tt, nanmean(all_high10piezo_preomit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_high10piezo_preomit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low10piezo_preomit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_low10piezo_preomit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice, ' num2str(totIC_area) ' dendrites- High10 (black) vs low10 (blue) pre-reward movements']);
        savefig([share_out '\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_HighVsLow10PreRewPiezo_' mouse_str '.fig'])     
        
        figure;
        subplot(2,2,1)
        shadedErrorBar(tt, nanmean(all_high25piezo_postrew(:,ind),2).*(1000./frameRateHz), (nanstd(all_high25piezo_postrew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low25piezo_postrew(:,ind),2).*(1000./frameRateHz), (nanstd(all_low25piezo_postrew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Rewarded trials')
        ylim([-1 12])
        vline(600)
        
        subplot(2,2,3)
        shadedErrorBar(tt, nanmean(all_high25piezo_postrew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_high25piezo_postrew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low25piezo_postrew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_low25piezo_postrew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        
        subplot(2,2,2)
        shadedErrorBar(tt, nanmean(all_high25piezo_postomit(:,ind),2).*(1000./frameRateHz), (nanstd(all_high25piezo_postomit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low25piezo_postomit(:,ind),2).*(1000./frameRateHz), (nanstd(all_low25piezo_postomit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omission trials')
        ylim([-1 12])
        
        subplot(2,2,4)
        shadedErrorBar(tt, nanmean(all_high25piezo_postomit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_high25piezo_postomit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low25piezo_postomit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_low25piezo_postomit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice, ' num2str(totIC_area) ' dendrites- High25 (black) vs low25 (blue) post-reward movements']);
        savefig([share_out '\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_HighVsLow25PostRewPiezo_' mouse_str '.fig'])     
    
        figure;
        subplot(2,2,1)
        shadedErrorBar(tt, nanmean(all_high10piezo_postrew(:,ind),2).*(1000./frameRateHz), (nanstd(all_high10piezo_postrew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low10piezo_postrew(:,ind),2).*(1000./frameRateHz), (nanstd(all_low10piezo_postrew(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Rewarded trials')
        ylim([-1 12])
        vline(600)
        
        subplot(2,2,3)
        shadedErrorBar(tt, nanmean(all_high10piezo_postrew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_high10piezo_postrew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low10piezo_postrew_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_low10piezo_postrew_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        
        subplot(2,2,2)
        shadedErrorBar(tt, nanmean(all_high10piezo_postomit(:,ind),2).*(1000./frameRateHz), (nanstd(all_high10piezo_postomit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low10piezo_postomit(:,ind),2).*(1000./frameRateHz), (nanstd(all_low10piezo_postomit(:,ind),[],2)./sqrt(totIC_area)).*(1000./frameRateHz),'b');
        xlabel('Time from cue')
        ylabel('Spike rate (Hz)')
        title('Omission trials')
        ylim([-1 12])
        
        subplot(2,2,4)
        shadedErrorBar(tt, nanmean(all_high10piezo_postomit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_high10piezo_postomit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'k');
        hold on
        shadedErrorBar(tt, nanmean(all_low10piezo_postomit_piezo(:,find(expt_areas(i,:))),2), (nanstd(all_low10piezo_postomit_piezo(:,find(expt_areas(i,:))),[],2)./sqrt(sum(expt_areas(i,:),2))),'b');
        hold on
        xlabel('Time from cue')
        ylabel('Piezo (V)')
        suptitle(['Day ' num2str(id) ' Area' area_list{i} ': n= ' num2str(totExp(i)) ' mice, ' num2str(totIC_area) ' dendrites- High10 (black) vs low10 (blue) post-reward movements']);
        savefig([share_out '\CC_summary\Day' num2str(id) '\Day' num2str(id) '_Summary_Area' area_list{i} '_HighVsLow10PostRewPiezo_' mouse_str '.fig'])     
    
    end
end
    
