% CRP_meta_k_means
%script for running meta k-means analysis on the 2-photon data from the CRP
%experiments.

%Since we are using the binary reports of event times we can just string
%together trials without worrying about the spliced ends of the trials due
%to laser power modulation.

clear all
close all
data_dir =  'Y:\home\jake\Analysis\Cue_reward_pairing_analysis\CC_analysis_2P';
for id = [1];
    
    %toggle for LS vs Crus datasets
    CRP_OT_expt_list_LS_jh;
    %CRP_OT_expt_list_Crus_jh;
    
    nexp = size(expt(id).date,1);
    fprintf(['Day: ' num2str(id) '\n'])
    
    for iexp = 1%1:nexp
        %get mouse ID info
        mouse = strtrim(expt(id).mouse(iexp,:));
        date = expt(id).date(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];
        
        %load data  _targetAlign
        load(fullfile(data_dir, img_fn, [img_fn '_targetAlign.mat']), 'targetAlign_events', 'targetAligndFoverF', 'frameRateHz', 'prewin_frames', 'postwin_frames');
        
        %if anaylzing df/f traces
        img_mat = targetAligndFoverF;
        %if analyzing calcium events
%        img_mat = targetAlign_events;
        
        img_mat = img_mat([3:end],:,[2:end]); %first trial is often all NaNs   %first two frames are often NaNs
        eventsmat = reshape(img_mat, size(img_mat,2), (size(img_mat,1)*size(img_mat,3))); %dim1=dend dim2=frames
        events_nan_ind = find(isnan( sum(eventsmat) )); %look across all frames to find any NaNs. 
        eventsmat(:, [events_nan_ind]) = []; %if a given frame has NaNs then remove it
        
        %perform meta k-means clustering 
        [allclusters, centroidcorr, dendmem, dunnsinitial] = meta_k_means_2P(eventsmat);
        sc = 5;
        [allclusters, centroidcorr, dendmem, dunnsinitial] = meta_k_means(eventsmat);
        %save and plot outputs
        
    end
end