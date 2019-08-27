% main2 
% run analysis for extracted tc
% 1. load behavior data
% 2. getTC_events extract tc based on trial 
% 3. TC_quantification find cell to each trial type and do within session
% analysis
% 4. getTC_Spike extract event info based on trial
% 5. spike_quantification finds firing rate and spikes
clear

file_info;
usFacs = 100;
%figure;
n_sub = zeros(30,6);
for sub = 1:30                         % [29 ,30]%[8, 11, 13, 15:20] %36
    for rID = 1:2
        file_info;
        mouse = mouseID{sub}
        dateID = date;
        date = dateID{sub}
        data_dir = fullfile('Z:\Analysis\2P Analysis\Lever');
        
        
        output_dest = fullfile('\\crash.dhe.duke.edu\data\home\jake\2P_Analysis\2P_Analysis',[date, '_', runID{rID}, '_', mouse],'\');
%         if ~exist(output_dest)
%             mkdir(output_dest);
%         end
        
        %subfold = fullfile(data_dir,[dateID{sub} '_' mouseID{sub} '_run' runID{rID} '\']);
        subfold = fullfile('\\crash.dhe.duke.edu\data\home\jake\2P_Analysis\2P_Analysis',[dateID{sub}, '_', runID{rID}, '_', mouseID{sub}],'\');
        dataName   = dir(fullfile(subfold,'*frame_times.mat'));
        
        behav_dir = '\\crash.dhe.duke.edu\data\home\andrew\Behavior\Data\';
        tc_dir  = fullfile('\\crash.dhe.duke.edu\data\home\jake\2P_Analysis\2P_Analysis',[dateID{sub}, '_', runID{rID}, '_', mouseID{sub}],'\');
       
        dest = tc_dir;
        dest_sub = dest;
        
        if exist(dest_sub)
            %load([tc_dir, 'ROI_TCs.mat']);
            
            dataName   = dir([behav_dir, '*', behavID{sub}, '-', date(1:6), '*']);
            load([behav_dir, dataName(end).name]);
            
            getTC_events;
            verify_trial_type_with_RTs
            TC_quantification;
            close all
            
            %subplot(5,6,sub)
            RTs = celleqel2mat_padded(input.reactTimesMs);
            RT_c = RTs;
            RT_c(find(RTs<200)) = nan;
            RT_c(find(RTs>950)) = nan;
            edges = 200:75:600;
            n_sub(sub,:) = histc(RT_c,edges);
            
            getTC_Spike
            spike_quantification
%             close all
        end
       
        clearvars -except mouseID sub rID n_sub
    end
end
