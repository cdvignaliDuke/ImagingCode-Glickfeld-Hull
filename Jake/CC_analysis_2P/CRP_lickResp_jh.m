clear all
close all
CRP_expt_list_all

lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';

for id = 1
    nexp = size(expt(id).date,1);
    for iexp = 1
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

    end
end

            
            
    