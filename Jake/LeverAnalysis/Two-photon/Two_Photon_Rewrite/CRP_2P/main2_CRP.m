%rewrite of main2_Cuepair to work with new datasets and analyses. Rewrite starting 8/19/18 JH

clear
file_info_CRP_all;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration';

for session_subset = 2
    %label session subsets to use
    if session_subset ==1
        days_subset = days_1;
        nPCA_to_use = nPCA_to_use_1;
    elseif session_subset ==2
        days_subset = days_post;
        nPCA_to_use = nPCA_to_use_post;
    elseif session_subset ==3
        days_subset = days_UR;
        nPCA_to_use = nPCA_to_use_UR;
    elseif session_subset ==4
        days_subset = days_1000;
        nPCA_to_use = nPCA_to_use_1000;
    end
    
    for sub = [1:length(days_subset)];
        if  isempty(days_subset{sub}) 
            continue
        end
        %collect session/mouse information
        session = days_subset{sub};
        session_date = days_subset{sub}(1:6);
        if session(end-2) == 'g'
            mouse_num = ['9', session(end-1:end)];
            mouse_ID = ['img', session(end-1:end)];
        elseif session(end-2) =='0'
            mouse_num = session(end-2:end);
            mouse_ID = ['img', session(end-2:end)];
        end
        session
        
        %set pathnames
        data_2P_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
        for rID = 1:2
            if exist([data_2P_dir, mouse_ID, '_000_', runID{rID}, '.sbx'], 'file') == 2
                break
            end
        end
        if strcmp(session, '170426_img91') | strcmp(session, '170522_img89')
            rID=2;
        end
        crp_main_dir =  fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID], '\');
        crp_spikes_dir = fullfile(crp_main_dir, 'spike_outputs\');
%         crp_spikes_dir = fullfile(crp_main_dir, 'spike_outputs2\');
        if ~exist(crp_spikes_dir)
            mkdir(crp_spikes_dir);
        end
        dest = crp_main_dir;
        dest_sub = crp_spikes_dir;
        
        %load behavior data
        dataName   = dir([behav_dir, '*', 'i', mouse_num, '-',session_date, '*']);
        load([behav_dir, dataName(end).name]);
        
        %getTC_Spike_CRP; %umcomment to run, extract spikes
        spike_quantification_CRP2; %umcomment to run, quantify spikes for rates and PSTH
        %plot_session_PSTH2;
    end
end

%Summary_events_1_CRP2;


