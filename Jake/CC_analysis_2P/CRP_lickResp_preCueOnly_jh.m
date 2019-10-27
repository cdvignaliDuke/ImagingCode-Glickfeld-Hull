clear all
close all
%CRP_expt_list_all
CRP_OT_expt_list_Crus_jh
lg_out = 'Y:\home\jake\Analysis\Cue_reward_pairing_analysis\CC_analysis_2P';

for id = 2:6
    nexp = size(expt(id).date,1);
    for iexp = 1:nexp
        mouse = strtrim(expt(id).mouse(iexp,:));
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];

        load(fullfile(lg_out,img_fn, [img_fn '_input.mat']))
        load(fullfile(lg_out,img_fn, [img_fn '_targetAlign.mat']))
        load(fullfile(lg_out,img_fn, [img_fn '_cueAlignLick.mat']))

        nIC = size(targetAlign_events,2);
        nTrials = size(targetAlign_events,3);
        singleLick_frames = round(0.25.*frameRateHz);
        tl = (-singleLick_frames:singleLick_frames-1).*(1000./frameRateHz);
        precue_single_lick = [];
        precue_single_lick_df = [];
        for itrial = 1:nTrials
            ind = find(lickCueAlign(1:prewin_frames,itrial));
            if length(ind) == 1 & ind>singleLick_frames & prewin_frames-ind>singleLick_frames
                precue_single_lick = cat(3,precue_single_lick, targetAlign_events(ind-singleLick_frames:ind+singleLick_frames-1,:,itrial));
                precue_single_lick_df = cat(3,precue_single_lick_df, targetAligndFoverF(ind-singleLick_frames:ind+singleLick_frames-1,:,itrial));
            elseif length(ind)>1
                for i = 1:length(ind)
                    if i == 1 
                        if ind(i)>singleLick_frames & prewin_frames-ind(i)>singleLick_frames & ind(i+1)-ind(i)>=singleLick_frames
                            precue_single_lick = cat(3,precue_single_lick, targetAlign_events(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            precue_single_lick_df = cat(3,precue_single_lick_df, targetAligndFoverF(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                        end
                    elseif i>1 & length(ind) == i 
                        if prewin_frames-ind(i)>singleLick_frames & ind(i)>singleLick_frames & ind(i)-ind(i-1)>=singleLick_frames
                            precue_single_lick = cat(3,precue_single_lick, targetAlign_events(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            precue_single_lick_df = cat(3,precue_single_lick_df, targetAligndFoverF(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                        end
                    elseif i>1 & length(ind) > i 
                        if ind(i)>singleLick_frames &  ind(i)-ind(i-1)>=singleLick_frames & ind(i+1)-ind(i)>=singleLick_frames & prewin_frames-ind(i)>singleLick_frames
                            precue_single_lick = cat(3,precue_single_lick, targetAlign_events(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            precue_single_lick_df = cat(3,precue_single_lick_df, targetAligndFoverF(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                        end
                    end
                end
            end
        end
        precue_lick_burst = [];
        precue_lick_burst_df = [];
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
                            precue_lick_burst = cat(3,precue_lick_burst, targetAlign_events(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            precue_lick_burst_df = cat(3,precue_lick_burst_df, targetAligndFoverF(ind(i)-singleLick_frames:ind(i)+singleLick_frames-1,:,itrial));
                            burst_trial = [burst_trial itrial];
                            burst_ind = [burst_ind ind(i)];
                        end
                    end
                end
            end
        end
        if size(precue_lick_burst,1)==0
            n_burst = 0;
        else
            n_burst = size(precue_lick_burst,3);
        end
        if size(precue_single_lick,1)==0
            n_lick = 0;
        else
            n_lick = size(precue_single_lick,3);
        end
        
        figure; 
        if n_lick>0
            subplot(2,2,1)
            shadedErrorBar(tl, mean(mean(precue_single_lick_df,2),3), std(mean(precue_single_lick_df,3),[],2)./sqrt(nIC));
            ylabel('dF/F')
            xlabel('Time from lick')
            ylim([-0.1 0.5])
            title(['Single lick'])
            subplot(2,2,2)
            shadedErrorBar(tl, mean(mean(precue_single_lick,2),3).*(1000./frameRateHz), std(mean(precue_single_lick,3),[],2)./sqrt(nIC).*(1000./frameRateHz));
            ylabel('Spike rate (Hz)')
            xlabel('Time from lick')
            ylim([0 10])
            title(['Single lick'])
        end
           
        if n_burst>0
            subplot(2,2,3)
            shadedErrorBar(tl, mean(mean(precue_lick_burst_df,2),3), std(mean(precue_lick_burst_df,3),[],2)./sqrt(nIC));
            ylabel('dF/F')
            xlabel('Time from lick')
            ylim([-0.1 0.5])
            title(['Lick burst'])
            subplot(2,2,4)
            shadedErrorBar(tl, mean(mean(precue_lick_burst,2),3).*(1000./frameRateHz), std(mean(precue_lick_burst,3),[],2)./sqrt(nIC).*(1000./frameRateHz));
            ylabel('Spike rate (Hz)')
            xlabel('Time from lick')
            ylim([0 10])
            title(['Lick burst'])
        end
            
        suptitle([mouse ' ' date ' pre-Cue lick response: single (' num2str(n_lick) '); burst (' num2str(n_burst) ')'])
        savefig(fullfile(lg_out,img_fn, [img_fn '_lickAlignSpiking_preCue_allCells.fig']))
        tl_precue = tl;
        save(fullfile(lg_out,img_fn, [img_fn '_lickResp_preCue.mat']), 'precue_single_lick', 'precue_lick_burst', 'precue_single_lick_df', 'precue_lick_burst_df', 'frameRateHz', 'tl_precue')

    end
end

            
            
    