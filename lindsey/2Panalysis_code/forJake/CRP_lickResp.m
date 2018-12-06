clear all
close all
CRP_expt_list
id = 1;
lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';
nexp = size(expt(id).date,1);
for iexp = 1:6
    mouse = expt(id).mouse(iexp,:);
    date = expt(id).date(iexp,:);
    run = expt(id).run(iexp,:);
    fprintf([date ' ' mouse '\n'])
    img_fn = [date '_' mouse];
    
    load(fullfile(lg_out,img_fn, [img_fn '_input.mat']))
    load(fullfile(lg_out,img_fn, [img_fn '_targetAlign.mat']))
    load(fullfile(lg_out,img_fn, [img_fn '_cueAlignLick.mat']))
   
    nIC = size(targetAlign_events,2);
    nTrials = size(targetAlign_events,3);
    singleLick_frames = round(0.5.*frameRateHz);
    tl = (-singleLick_frames:singleLick_frames-1).*(1000./frameRateHz);
    single_lick = [];
    single_lick_df = [];
    for itrial = 1:nTrials
        ind = find(lickCueAlign(1:prewin_frames,itrial));
        if length(ind) == 1 & ind>singleLick_frames & prewin_frames-ind>singleLick_frames
        	single_lick = cat(3,single_lick, targetAlign_events(ind-singleLick_frames:ind+singleLick_frames-1,:,itrial));
            single_lick_df = cat(3,single_lick_df, targetAligndFoverF(ind-singleLick_frames:ind+singleLick_frames-1,:,itrial));
        elseif length(ind)>1
            for i = 1:length(ind)
                if i == 1 
                    if ind(i)>singleLick_frames & prewin_frames-ind(i)>singleLick_frames & ind(i+1)-ind(i)>=singleLick_frames
                        single_lick = cat(3,single_lick, targetAlign_events(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                        single_lick_df = cat(3,single_lick_df, targetAligndFoverF(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                    end
                elseif i>1 & length(ind) == i 
                    if prewin_frames-ind(i)>singleLick_frames & ind(i)>singleLick_frames & ind(i)-ind(i-1)>=singleLick_frames
                        single_lick = cat(3,single_lick, targetAlign_events(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                        single_lick_df = cat(3,single_lick_df, targetAligndFoverF(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                    end
                elseif i>1 & length(ind) > i 
                    if ind(i)>singleLick_frames &  ind(i)-ind(i-1)>=singleLick_frames & ind(i+1)-ind(i)>=singleLick_frames & prewin_frames-ind(i)>singleLick_frames
                        single_lick = cat(3,single_lick, targetAlign_events(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                        single_lick_df = cat(3,single_lick_df, targetAligndFoverF(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                    end
                end
            end
        end
    end
    lick_burst = [];
    lick_burst_df = [];
    burst_trial = [];
    burst_ind = [];
    for itrial = 1:nTrials
        ind = find(lickCueAlign(1:prewin_frames,itrial));
        if length(ind)>=3 & ind(1)>=singleLick_frames
            for i = 1:length(ind)-2    
                if ind(i+2)-ind(i)<= lickSearch_frames & prewin_frames-ind(i)>singleLick_frames & ind(i)>singleLick_frames
                    if find(burst_trial == itrial) & ind(i)-burst_ind(end)<singleLick_frames
                        continue
                    else
                        lick_burst = cat(3,lick_burst, targetAlign_events(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                        lick_burst_df = cat(3,lick_burst_df, targetAligndFoverF(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                        burst_trial = [burst_trial itrial];
                        burst_ind = [burst_ind ind(i)];
                    end
                end
            end
        end
    end
    if size(lick_burst,1)==0
        n_burst = 0;
    else
        n_burst = size(lick_burst,3);
    end
    if size(single_lick,1)==0
        n_lick = 0;
    else
        n_lick = size(single_lick,3);
    end
    h_burst = zeros(1,nIC);
    h_lick = zeros(1,nIC);
    p_burst = zeros(1,nIC);
    p_lick = zeros(1,nIC);
    h_burst_df = zeros(1,nIC);
    h_lick_df = zeros(1,nIC);
    p_burst_df = zeros(1,nIC);
    p_lick_df = zeros(1,nIC);
    for iC = 1:nIC
        if n_burst>0
            [h_burst(:,iC) p_burst(:,iC)] = ttest(squeeze(mean(lick_burst(6:10,iC,:),1)),squeeze(mean(lick_burst(11:15,iC,:),1)),'tail','left'); 
            [h_burst_df(:,iC) p_burst_df(:,iC)] = ttest(squeeze(mean(lick_burst_df(6:10,iC,:),1)),squeeze(mean(lick_burst_df(11:15,iC,:),1)),'tail','left'); 
        end
        if n_lick>0
            [h_lick(:,iC) p_lick(:,iC)] = ttest(squeeze(mean(single_lick(6:10,iC,:),1)),squeeze(mean(single_lick(11:15,iC,:),1)),'tail','left'); 
            [h_lick_df(:,iC) p_lick_df(:,iC)] = ttest(squeeze(mean(single_lick_df(6:10,iC,:),1)),squeeze(mean(single_lick_df(11:15,iC,:),1)),'tail','left'); 
        end
    end
    
    %corrects for cases with no spikes in either window
    h_burst(isnan(h_burst)) = 0;
    h_lick(isnan(h_lick)) = 0;
    h_burst_df(isnan(h_burst_df)) = 0;
    h_lick_df(isnan(h_lick_df)) = 0;
    
    figure; 
    if n_lick>0
        subplot(2,3,1)
        shadedErrorBar(tl, mean(mean(single_lick,2),3).*(1000./frameRateHz), std(mean(single_lick,3),[],2)./sqrt(nIC).*(1000./frameRateHz));
        ylabel('Spike rate (Hz)')
        xlabel('Time from lick')
        ylim([0 10])
        title(['Single lick- n= ' num2str(sum(h_lick,2)) '/' num2str(nIC) ' resp'])
        subplot(2,3,2)
        shadedErrorBar(tl, mean(mean(single_lick(:,find(h_lick),:),2),3).*(1000./frameRateHz), std(mean(single_lick(:,find(h_lick),:),3),[],2)./sqrt(nIC).*(1000./frameRateHz));
        ylabel('Spike rate (Hz)')
        xlabel('Time from lick')
        ylim([0 10])
        title(['Single lick resp'])
        subplot(2,3,3)
        shadedErrorBar(tl, mean(mean(single_lick(:,find(~h_lick),:),2),3).*(1000./frameRateHz), std(mean(single_lick(:,find(~h_lick),:),3),[],2)./sqrt(nIC).*(1000./frameRateHz));
        ylabel('Spike rate (Hz)')
        xlabel('Time from lick')
        ylim([0 10])
        title(['Single lick not resp'])
    end
    if n_burst>0
        subplot(2,3,4)
        shadedErrorBar(tl, mean(mean(lick_burst,2),3).*(1000./frameRateHz), std(mean(lick_burst,3),[],2)./sqrt(nIC).*(1000./frameRateHz));
        ylabel('Spike rate (Hz)')
        xlabel('Time from lick')
        ylim([0 10])
        title(['Lick burst- n= ' num2str(sum(h_burst,2)) '/' num2str(nIC) ' resp'])
        subplot(2,3,5)
        shadedErrorBar(tl, mean(mean(lick_burst(:,find(h_burst),:),2),3).*(1000./frameRateHz), std(mean(lick_burst(:,find(h_burst),:),3),[],2)./sqrt(nIC).*(1000./frameRateHz));
        ylabel('Spike rate (Hz)')
        xlabel('Time from lick')
        ylim([0 10])
        title(['Lick burst resp'])
        subplot(2,3,6)
        shadedErrorBar(tl, mean(mean(lick_burst(:,find(~h_burst),:),2),3).*(1000./frameRateHz), std(mean(lick_burst(:,find(~h_burst),:),3),[],2)./sqrt(nIC).*(1000./frameRateHz));
        ylabel('Spike rate (Hz)')
        xlabel('Time from lick')
        ylim([0 10])
        title(['Lick burst not resp'])
    end
    suptitle([mouse ' ' date ' pre-Cue lick response: single (' num2str(n_lick) '); burst (' num2str(n_burst) ')'])
    savefig(fullfile(lg_out,img_fn, [img_fn '_lickAlignSpiking_preCue.fig']))
    
    figure; 
    if n_lick>0
        subplot(2,3,1)
        shadedErrorBar(tl, mean(mean(single_lick_df,2),3), std(mean(single_lick_df,3),[],2)./sqrt(nIC));
        ylabel('dF/F')
        xlabel('Time from lick')
        ylim([-0.1 0.5])
        title(['Single lick- n= ' num2str(sum(h_lick_df,2)) '/' num2str(nIC) ' resp'])
        subplot(2,3,2)
        shadedErrorBar(tl, mean(mean(single_lick_df(:,find(h_lick_df),:),2),3), std(mean(single_lick_df(:,find(h_lick_df),:),3),[],2)./sqrt(nIC));
        ylabel('dF/F')
        xlabel('Time from lick')
        ylim([-0.1 0.5])
        title(['Single lick resp'])
        subplot(2,3,3)
        shadedErrorBar(tl, mean(mean(single_lick_df(:,find(~h_lick_df),:),2),3), std(mean(single_lick_df(:,find(~h_lick_df),:),3),[],2)./sqrt(nIC));
        ylabel('dF/F')
        xlabel('Time from lick')
        ylim([-0.1 0.5])
        title(['Single lick not resp'])
    end
    if n_burst>0
        subplot(2,3,4)
        shadedErrorBar(tl, mean(mean(lick_burst_df,2),3), std(mean(lick_burst_df,3),[],2)./sqrt(nIC));
        ylabel('dF/F')
        xlabel('Time from lick')
        ylim([-0.1 0.5])
        title(['Lick burst- n= ' num2str(sum(h_burst_df,2)) '/' num2str(nIC) ' resp'])
        subplot(2,3,5)
        shadedErrorBar(tl, mean(mean(lick_burst_df(:,find(h_burst_df),:),2),3), std(mean(lick_burst_df(:,find(h_burst_df),:),3),[],2)./sqrt(nIC));
        ylabel('dF/F')
        xlabel('Time from lick')
        ylim([-0.1 0.5])
        title(['Lick burst resp'])
        subplot(2,3,6)
        shadedErrorBar(tl, mean(mean(lick_burst_df(:,find(~h_burst_df),:),2),3), std(mean(lick_burst_df(:,find(~h_burst_df),:),3),[],2)./sqrt(nIC));
        ylabel('dF/F')
        xlabel('Time from lick')
        ylim([-0.1 0.5])
        title(['Lick burst not resp'])
    end
    suptitle([mouse ' ' date ' pre-Cue lick response: single (' num2str(n_lick) '); burst (' num2str(n_burst) ')'])
    savefig(fullfile(lg_out,img_fn, [img_fn '_lickAlignDFoverF_preCue.fig']))
    
    save(fullfile(lg_out,img_fn, [img_fn '_lickResp.mat']), 'single_lick', 'lick_burst', 'single_lick_df', 'lick_burst_df', 'singleLick_frames', 'frameRateHz', 'tl', 'h_burst', 'h_lick','h_burst_df', 'h_lick_df')

    if sum(ismember(1:6, iexp))

        load(fullfile(lg_out,img_fn, [img_fn '_ROI_TCs.mat']))
        load(fullfile(lg_out,img_fn, [img_fn '_ROI_spikes.mat']))
        
        cTargetOn = celleqel2mat_padded(input.cTargetOn);

        pretarget_frames_1 = round(15000./frameRateHz);
        pretarget_frames_2 = round(0./frameRateHz);
        itiAlign_tc = nan(pretarget_frames_1-pretarget_frames_2,nIC,nTrials);
        itiAlign_events = nan(pretarget_frames_1-pretarget_frames_2,nIC,nTrials);
        totWindow = pretarget_frames_1-pretarget_frames_2;
        itiLicks = nan(totWindow,nTrials);
        for itrial = 1:nTrials
            counterVals = input.counterValues{itrial};
            itiAlign_tc(:,:,itrial) = tc_avg(cTargetOn(itrial)-pretarget_frames_1:cTargetOn(itrial)-pretarget_frames_2-1,:);
            itiAlign_events(:,:,itrial) = all_events(cTargetOn(itrial)-pretarget_frames_1:cTargetOn(itrial)-pretarget_frames_2-1,:);
            itiLicks(:,itrial) =lickTC{itrial}(1,cTargetOn(itrial)-pretarget_frames_1-counterVals(1):cTargetOn(itrial)+pretarget_frames_2-1-counterVals(1))';
        end
        itiAlign_dFoverF = (itiAlign_tc-mean(itiAlign_tc,1))./mean(itiAlign_tc,1);

        singleLick_frames = round(0.5.*frameRateHz);
        tl = (-singleLick_frames:singleLick_frames-1).*(1000./frameRateHz);
        iti_single_lick = [];
        iti_single_lick_df = [];
        iti_lick_ind = [];
        for itrial = 1:nTrials
            ind = find(itiLicks(:,itrial));
            if length(ind) == 1 & ind>singleLick_frames & totWindow-ind>singleLick_frames
                iti_single_lick = cat(3,iti_single_lick, itiAlign_events(ind-singleLick_frames:ind+singleLick_frames-1,:,itrial));
                iti_single_lick_df = cat(3,iti_single_lick_df, itiAlign_dFoverF(ind-singleLick_frames:ind+singleLick_frames-1,:,itrial));
                iti_lick_ind = [iti_lick_ind ind];
            elseif length(ind)>1
                for i = 1:length(ind)
                    if i == 1 
                        if ind(i)>singleLick_frames & totWindow-ind(i)>singleLick_frames & ind(i+1)-ind(i)>=singleLick_frames
                            iti_single_lick = cat(3,iti_single_lick, itiAlign_events(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            iti_single_lick_df = cat(3,iti_single_lick_df, itiAlign_dFoverF(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            iti_lick_ind = [iti_lick_ind ind(i)];
                        end
                    elseif i>1 & length(ind) == i 
                        if totWindow-ind(i)>singleLick_frames & ind(i)>singleLick_frames & ind(i)-ind(i-1)>=singleLick_frames
                            iti_single_lick = cat(3,iti_single_lick, itiAlign_events(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            iti_single_lick_df = cat(3,iti_single_lick_df, itiAlign_dFoverF(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            iti_lick_ind = [iti_lick_ind ind(i)];
                        end
                    elseif i>1 & length(ind) > i 
                        if ind(i)>singleLick_frames &  ind(i)-ind(i-1)>=singleLick_frames & ind(i+1)-ind(i)>=singleLick_frames & totWindow-ind(i)>singleLick_frames
                            iti_single_lick = cat(3,iti_single_lick, itiAlign_events(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            iti_single_lick_df = cat(3,iti_single_lick_df, itiAlign_dFoverF(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            iti_lick_ind = [iti_lick_ind ind(i)];
                        end
                    end
                end
            end
        end
        iti_lick_burst = [];
        iti_lick_burst_df = [];
        iti_burst_trial = [];
        iti_burst_ind = [];
        for itrial = 1:nTrials
            ind = find(itiLicks(:,itrial));
            if length(ind)>=3 & ind(1)>=singleLick_frames
                for i = 1:length(ind)-2    
                    if ind(i+2)-ind(i)<= lickSearch_frames & totWindow-ind(i)>singleLick_frames & ind(i)>singleLick_frames
                        if find(iti_burst_trial == itrial) & ind(i)-iti_burst_ind(end)<singleLick_frames
                            continue
                        else
                            iti_lick_burst = cat(3,iti_lick_burst, itiAlign_events(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            iti_lick_burst_df = cat(3,iti_lick_burst_df, itiAlign_dFoverF(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            iti_burst_trial = [iti_burst_trial itrial];
                            iti_burst_ind = [iti_burst_ind ind(i)];
                        end
                    end
                end
            end
        end
        if size(iti_lick_burst,1)==0
            iti_n_burst = 0;
        else
            iti_n_burst = size(iti_lick_burst,3);
        end
        if size(iti_single_lick,1)==0
            iti_n_lick = 0;
        else
            iti_n_lick = size(iti_single_lick,3);
        end
        h_itiburst = zeros(1,nIC);
        h_itilick = zeros(1,nIC);
        p_itiburst = zeros(1,nIC);
        p_itilick = zeros(1,nIC);
        h_itiburst_df = zeros(1,nIC);
        h_itilick_df = zeros(1,nIC);
        p_itiburst_df = zeros(1,nIC);
        p_itilick_df = zeros(1,nIC);
        for iC = 1:nIC
            if iti_n_burst>3
                [h_itiburst(:,iC) p_itiburst(:,iC)] = ttest(squeeze(mean(iti_lick_burst(6:10,iC,:),1)),squeeze(mean(iti_lick_burst(11:15,iC,:),1)),'tail','left'); 
                [h_itiburst_df(:,iC) p_itiburst_df(:,iC)] = ttest(squeeze(mean(iti_lick_burst_df(6:10,iC,:),1)),squeeze(mean(iti_lick_burst_df(11:15,iC,:),1)),'tail','left'); 
            end
            if n_lick>3
                [h_itilick(:,iC) p_itilick(:,iC)] = ttest(squeeze(mean(iti_single_lick(6:10,iC,:),1)),squeeze(mean(iti_single_lick(11:15,iC,:),1)),'tail','left'); 
                [h_itilick_df(:,iC) p_itilick_df(:,iC)] = ttest(squeeze(mean(iti_single_lick_df(6:10,iC,:),1)),squeeze(mean(iti_single_lick_df(11:15,iC,:),1)),'tail','left'); 
            end
        end

        %corrects for cases with no spikes in either window
        h_itiburst(isnan(h_itiburst)) = 0;
        h_itilick(isnan(h_itilick)) = 0;
        h_itiburst_df(isnan(h_itiburst_df)) = 0;
        h_itilick_df(isnan(h_itilick_df)) = 0;

        figure; 
        if iti_n_lick>0
            subplot(2,3,1)
            shadedErrorBar(tl, mean(mean(iti_single_lick,2),3).*(1000./frameRateHz), std(mean(iti_single_lick,3),[],2)./sqrt(nIC).*(1000./frameRateHz));
            ylabel('Spike rate (Hz)')
            xlabel('Time from lick')
            ylim([0 10])
            title(['Single lick- n= ' num2str(sum(h_itilick,2)) '/' num2str(nIC) ' resp'])
            subplot(2,3,2)
            shadedErrorBar(tl, mean(mean(iti_single_lick(:,find(h_itilick),:),2),3).*(1000./frameRateHz), std(mean(iti_single_lick(:,find(h_itilick),:),3),[],2)./sqrt(nIC).*(1000./frameRateHz));
            ylabel('Spike rate (Hz)')
            xlabel('Time from lick')
            ylim([0 10])
            title(['Single lick resp'])
            subplot(2,3,3)
            shadedErrorBar(tl, mean(mean(iti_single_lick(:,find(~h_itilick),:),2),3).*(1000./frameRateHz), std(mean(iti_single_lick(:,find(~h_itilick),:),3),[],2)./sqrt(nIC).*(1000./frameRateHz));
            ylabel('Spike rate (Hz)')
            xlabel('Time from lick')
            ylim([0 10])
            title(['Single lick not resp'])
        end
        if iti_n_burst>0
            subplot(2,3,4)
            shadedErrorBar(tl, mean(mean(iti_lick_burst,2),3).*(1000./frameRateHz), std(mean(iti_lick_burst,3),[],2)./sqrt(nIC).*(1000./frameRateHz));
            ylabel('Spike rate (Hz)')
            xlabel('Time from lick')
            ylim([0 10])
            title(['Lick burst- n= ' num2str(sum(h_itiburst,2)) '/' num2str(nIC) ' resp'])
            subplot(2,3,5)
            shadedErrorBar(tl, mean(mean(iti_lick_burst(:,find(h_itiburst),:),2),3).*(1000./frameRateHz), std(mean(iti_lick_burst(:,find(h_itiburst),:),3),[],2)./sqrt(nIC).*(1000./frameRateHz));
            ylabel('Spike rate (Hz)')
            xlabel('Time from lick')
            ylim([0 10])
            title(['Lick burst resp'])
            subplot(2,3,6)
            shadedErrorBar(tl, mean(mean(iti_lick_burst(:,find(~h_itiburst),:),2),3).*(1000./frameRateHz), std(mean(iti_lick_burst(:,find(~h_itiburst),:),3),[],2)./sqrt(nIC).*(1000./frameRateHz));
            ylabel('Spike rate (Hz)')
            xlabel('Time from lick')
            ylim([0 10])
            title(['Lick burst not resp'])
        end
        suptitle([mouse ' ' date ' ITI lick response: single (' num2str(iti_n_lick) '); burst (' num2str(iti_n_burst) ')'])
        savefig(fullfile(lg_out,img_fn, [img_fn '_lickAlignSpiking_ITI.fig']))

        figure; 
        if iti_n_lick>0
            subplot(2,3,1)
            shadedErrorBar(tl, mean(mean(iti_single_lick_df,2),3), std(mean(iti_single_lick_df,3),[],2)./sqrt(nIC));
            ylabel('dF/F')
            xlabel('Time from lick')
            ylim([-0.1 0.5])
            title(['Single lick- n= ' num2str(sum(h_itilick_df,2)) '/' num2str(nIC) ' resp'])
            subplot(2,3,2)
            shadedErrorBar(tl, mean(mean(iti_single_lick_df(:,find(h_itilick_df),:),2),3), std(mean(iti_single_lick_df(:,find(h_itilick_df),:),3),[],2)./sqrt(nIC));
            ylabel('dF/F')
            xlabel('Time from lick')
            ylim([-0.1 0.5])
            title(['Single lick resp'])
            subplot(2,3,3)
            shadedErrorBar(tl, mean(mean(iti_single_lick_df(:,find(~h_itilick_df),:),2),3), std(mean(iti_single_lick_df(:,find(~h_itilick_df),:),3),[],2)./sqrt(nIC));
            ylabel('dF/F')
            xlabel('Time from lick')
            ylim([-0.1 0.5])
            title(['Single lick not resp'])
        end
        if iti_n_burst>0
            subplot(2,3,4)
            shadedErrorBar(tl, mean(mean(iti_lick_burst_df,2),3), std(mean(iti_lick_burst_df,3),[],2)./sqrt(nIC));
            ylabel('dF/F')
            xlabel('Time from lick')
            ylim([-0.1 0.5])
            title(['Lick burst- n= ' num2str(sum(h_itiburst_df,2)) '/' num2str(nIC) ' resp'])
            subplot(2,3,5)
            shadedErrorBar(tl, mean(mean(iti_lick_burst_df(:,find(h_itiburst_df),:),2),3), std(mean(iti_lick_burst_df(:,find(h_itiburst_df),:),3),[],2)./sqrt(nIC));
            ylabel('dF/F')
            xlabel('Time from lick')
            ylim([-0.1 0.5])
            title(['Lick burst resp'])
            subplot(2,3,6)
            shadedErrorBar(tl, mean(mean(iti_lick_burst_df(:,find(~h_itiburst_df),:),2),3), std(mean(iti_lick_burst_df(:,find(~h_itiburst_df),:),3),[],2)./sqrt(nIC));
            ylabel('dF/F')
            xlabel('Time from lick')
            ylim([-0.1 0.5])
            title(['Lick burst not resp'])
        end
        suptitle([mouse ' ' date ' ITI lick response: single (' num2str(iti_n_lick) '); burst (' num2str(iti_n_burst) ')'])
        savefig(fullfile(lg_out,img_fn, [img_fn '_lickAlign_dFoverF_ITI.fig']))

        edges = [1:ceil(totWindow/4):totWindow totWindow];
        [itilick_n itilick_ind] = histc(iti_lick_ind,edges);
        [itiburst_n itiburst_ind] = histc(iti_burst_ind,edges);
        figure;
        for i = 1:length(itilick_n)
            subplot(2,2,1)
            ind = find(itilick_ind == i);
            plot(tl, mean(mean(iti_single_lick(:,:,ind),2),3).*(1000./frameRateHz))
            ylabel('Spike Rate (Hz)')
            xlabel('Time from lick')
            title('Single lick')
            ylim([0 8])
            hold on;
            subplot(2,2,2)
            ind = find(itiburst_ind == i);
            plot(tl, mean(mean(iti_lick_burst(:,:,ind),2),3).*(1000./frameRateHz))
            ylabel('Spike Rate (Hz)')
            xlabel('Time from lick')
            title('Burst')
            ylim([0 8])
            hold on;
            subplot(2,2,3)
            ind = find(itilick_ind == i);
            plot(tl, mean(mean(iti_single_lick_df(:,:,ind),2),3).*(1000./frameRateHz))
            ylabel('dF/F')
            xlabel('Time from lick')
            title('Single lick')
            ylim([-2 8])
            hold on;
            subplot(2,2,4)
            ind = find(itiburst_ind == i);
            plot(tl, mean(mean(iti_lick_burst_df(:,:,ind),2),3).*(1000./frameRateHz))
            ylabel('dF/F')
            xlabel('Time from lick')
            title('Burst')
            ylim([-2 8])
            hold on;
        end
        suptitle([mouse ' ' date ' ITI lick response: single (' num2str(iti_n_lick) '); burst (' num2str(iti_n_burst) ')'])
        savefig(fullfile(lg_out,img_fn, [img_fn '_lickAlign_ITI_byTime.fig']))
        
        save(fullfile(lg_out,img_fn, [img_fn '_lickResp_ITI.mat']), 'iti_single_lick', 'iti_lick_burst', 'iti_single_lick_df', 'iti_lick_burst_df', 'singleLick_frames', 'frameRateHz', 'tl', 'h_itiburst', 'h_itilick','h_itiburst_df', 'h_itilick_df')
    end
    close all
end


            
            
    