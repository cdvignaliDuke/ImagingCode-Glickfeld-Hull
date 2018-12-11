% finds and maps responsive cells for each condition
% 1. calculate average timecourses for press/release events
% 2. calculate variability by trial over base and resp windows
% 3. ttest for significant responses
% 4. define cells by response to event

load([dest 'ROI_TCs.mat']);
load([dest_sub '_release_movies.mat'])
load([dest_sub '_press_movies.mat'])
load([dest_sub 'parse_behavior_LG.mat'])
load([dest_sub '_lick_stats_LG.mat']);
if ~isempty(lick_data)
    load([dest_sub '_lick_movies_LG.mat'])
end
%load([dest_sub 'parse_behavior.mat'])
nCells = size(success_movie,2);
n = ceil(sqrt(nCells));
if nCells <((n.^2)-n) %function
    n2= n-1;
else
    n2 = n;
end

%% 1. calculate average timecourses for press/release events
%average and sem across all release
% release_movie = cat(1,success_movie, fail_movie);
release_movie = release_long_movie;
avg_release = squeeze(mean(release_movie,1));
% avg_release_long = squeeze(mean(release_movie_long,1));
sem_release = squeeze(std(release_movie,1)./sqrt(size(release_movie,1)));
%average and sem across cells (and trials)
avg_release_all = mean(avg_release,1);
sem_release_all = std(avg_release,1)./sqrt(size(avg_release,1));

%average and sem across success
avg_success = squeeze(mean(success_movie,1));
avg_success = bsxfun(@minus, avg_success, avg_success(:,1));
avg_success_long = squeeze(mean(success_movie_long,1));
sem_success = squeeze(std(success_movie,1)./sqrt(size(avg_success,1)));
%average and sem across cells (and trials)
avg_success_all = mean(avg_success,1);
sem_success_all = std(avg_success,1)./sqrt(size(avg_success,1));
sem_success_long = squeeze(std(success_movie_long,1)./sqrt(size(avg_success_long,1)));

avg_early_success = squeeze(mean(early_success_movie,1));
sem_early_success = squeeze(std(early_success_movie,1)./sqrt(size(avg_early_success,1)));
avg_late_success = squeeze(mean(late_success_movie,1));
sem_late_success = squeeze(std(late_success_movie,1)./sqrt(size(avg_late_success,1)));

%average and sem across tooFast
avg_tooFast = squeeze(mean(tooFastCorrect_movie,1));
avg_tooFast = bsxfun(@minus, avg_tooFast, avg_tooFast(:,1));
avg_tooFast_long = squeeze(mean(tooFastCorrect_movie_long,1));
sem_tooFast = squeeze(std(tooFastCorrect_movie,1)./sqrt(size(tooFastCorrect_movie,1)));
%average and sem across cells (and trials)
avg_tooFast_all = mean(avg_tooFast,1);
sem_tooFast_all = std(avg_tooFast,1)./sqrt(size(avg_tooFast,1));

avg_early_tooFast = squeeze(mean(early_tooFast_movie,2));
avg_early_tooFast = bsxfun(@minus, avg_early_tooFast, avg_early_tooFast(:,1));
sem_early_tooFast = squeeze(std(early_tooFast_movie,1)./sqrt(size(avg_early_tooFast,1)));
avg_late_tooFast = squeeze(mean(late_tooFast_movie,2));
avg_late_tooFast = bsxfun(@minus, avg_late_tooFast, avg_late_tooFast(:,1));
sem_late_tooFast = squeeze(std(late_tooFast_movie,1)./sqrt(size(avg_late_tooFast,1)));

%average and sem across failures
avg_fail = squeeze(mean(fail_movie,1));
avg_fail = bsxfun(@minus, avg_fail, avg_fail(:,1));
avg_fail_long = squeeze(mean(fail_movie_long,1));
sem_fail = squeeze(std(fail_movie,1)./sqrt(size(avg_fail,1)));
%average and sem across ROIs (and trials)
avg_fail_all = mean(avg_fail,1);
sem_fail_all = std(avg_fail,1)./sqrt(size(avg_fail,1));

avg_early_fail = squeeze(nanmean(early_fail_movie,1));
avg_early_fail = bsxfun(@minus, avg_early_fail, avg_early_fail(:,1));
sem_early_fail = squeeze(nanstd(early_fail_movie,1)./sqrt(size(avg_early_fail,1)));

avg_late_fail = squeeze(nanmean(late_fail_movie,1));
avg_late_fail = bsxfun(@minus, avg_late_fail, avg_late_fail(:,1));
sem_late_fail = squeeze(nanstd(late_fail_movie,1)./sqrt(size(avg_late_fail,1)));


if ~isempty(lick_data)
    %average and sem across lickbout during ITI
    avg_lick = squeeze(mean(lick_movie,1));
    sem_lick = squeeze(std(lick_movie,1)./sqrt(size(avg_lick,1)));
    %average and sem across ROIs (and trials)
    avg_lick_all = mean(avg_lick,1);
    sem_lick_all = std(avg_lick,1)./sqrt(size(avg_lick,1));
    
    
    %average and sem across single lick during ITI
    avg_single_lick = squeeze(mean(single_lick_movie,1));
    sem_single_lick = squeeze(std(single_lick_movie,1)./sqrt(size(avg_single_lick,1)));
    %average and sem across ROIs (and trials)
    avg_single_lick_all = mean(avg_single_lick,1);
    sem_single_lick_all = std(avg_single_lick,1)./sqrt(size(avg_single_lick,1));

else
    avg_lick = [];
    avg_single_lick = [];
    sem_lick = []; sem_single_lick = [];
end

if ~isempty(omitReward_movie)
    avg_OR = squeeze(mean(omitReward_movie,1));
    avg_OR_long = squeeze(mean(omitReward_movie_long,1));
    sem_OR = squeeze(std(omitReward_movie,1)./sqrt(size(avg_OR,1)));
    %average and sem across ROIs (and trials)
    avg_OR_all = mean(avg_OR,1);
    sem_OR_all = std(avg_OR,1)./sqrt(size(avg_OR,1));
    sem_OR_long = squeeze(std(omitReward_movie_long,1)./sqrt(size(avg_OR_long,1)));
else
    avg_OR = []; avg_OR_long =[];
end

if ~isempty(itiReward_movie)
    avg_IR = squeeze(mean(itiReward_movie,1));
    avg_IR_long = squeeze(mean(itiReward_movie_long,1));
    sem_IR = squeeze(std(itiReward_movie,1)./sqrt(size(avg_IR,1)));
    %average and sem across ROIs (and trials)
    avg_IR_all = mean(avg_IR,1);
    sem_IR_all = std(avg_IR,1)./sqrt(size(avg_IR,1));
else
    avg_IR = []; avg_IR_long = [];
end

%average and sem across presses
avg_press = squeeze(mean(press_long_movie,1));
avg_press_long = squeeze(mean(press_long_movie_long,1));
sem_press = squeeze(std(press_long_movie,1)./sqrt(size(press_long_movie,1)));
%average and sem across ROIs (and trials)
avg_press_all = mean(avg_press,1);
sem_press_all = std(avg_press,1)./sqrt(size(avg_press,1));

save([dest_sub '_cell_TCs_LG.mat'], 'avg_release', 'avg_press','avg_success', 'avg_tooFast', 'avg_OR', ...
    'avg_fail', 'avg_press_long','avg_success_long', 'avg_tooFast_long', 'avg_fail_long', 'avg_OR_long', ...
    'sem_release', 'sem_press','sem_success', 'sem_tooFast', 'sem_fail', 'avg_lick', 'avg_single_lick', ...
    'sem_lick', 'sem_single_lick', 'avg_IR', 'avg_IR_long', 'avg_early_fail', 'avg_late_fail', ...
    'avg_early_success', 'avg_late_success', 'avg_early_tooFast', 'avg_late_tooFast')

% find rise time of release
peakIdx_release = findRisetime(cuerelease_movie, pre_cuerelease_frames)*double(ifi);
lever.cue2release(cue_remove_idx) = [];
peakIdx_cue = bsxfun(@plus, peakIdx_release, double(lever.cue2release'));

peak2release_std1 = std(abs(peakIdx_release),1); %all events individual cell 
peak2cue_std1 = std(peakIdx_cue,1);
peak2release_std2 = std(peakIdx_release, 0, 2); % all cells for each event
peak2cue_std2 = std(peakIdx_cue, 0, 2);

% save([dest_sub '_spike_variance.mat'], 'peak2release_std1', 'peak2cue_std1', 'peak2release_std2', 'peak2cue_std2');

%% 2. calculate response amplitude and variablity by trial over base and resp windows
%by release
base_release_window = 1:pre_release_frames-3; %selects a baseline window. Should rewrite this to use baseline times
base_release_window_new = pre_release_frames-2:pre_release_frames;
resp_release_window = pre_release_frames+1:pre_release_frames+round(100./double(ifi));  %selecting a response window in which to analyze the TC. Selects three consecutive frames where the lever event occurs during the first fraem

poss_succ_wins = [];  %using a sliding window to determine peak response to lever release 
poss_succ_win_vals = []; 
for i = pre_release_frames+1:pre_release_frames+1+round(300./double(ifi))
    poss_succ_win_vals = [poss_succ_win_vals; mean(mean(mean(success_movie(:,:,i-1:i+1)),3))];
    poss_succ_wins = [poss_succ_wins, i];
end
resp_success_window = poss_succ_wins(find(poss_succ_win_vals == max(poss_succ_win_vals))); %selects window with the peak response 
resp_success_window = [resp_success_window-1, resp_success_window, resp_success_window+1]; 

poss_fail_wins = [];  %using a sliding window to determine peak response to lever release 
poss_fail_win_vals = []; 
for i = pre_release_frames+1:pre_release_frames+1+round(300./double(ifi))
    poss_fail_win_vals = [poss_fail_win_vals; mean(mean(mean(fail_movie(:,:,i-1:i+1)),3))];
    poss_fail_wins = [poss_fail_wins, i];
end
fail_temp = (poss_fail_win_vals(poss_fail_win_vals>min(poss_fail_win_vals)));
% fail_vals = min(fail_temp(fail_temp>min(fail_temp)));
fail_vals = min(fail_temp);
resp_fail_window = poss_fail_wins(find(poss_fail_win_vals == fail_vals)); %selects window with the peak response 
resp_fail_window = [resp_fail_window-1, resp_fail_window, resp_fail_window+1]; 

poss_early_fail_wins = [];  %using a sliding window to determine peak response to lever release 
poss_early_fail_win_vals = []; 
for i = pre_release_frames+1:pre_release_frames+1+round(300./double(ifi))
    poss_early_fail_win_vals = [poss_early_fail_win_vals; mean(mean(mean(fail_movie(:,:,i-1:i+1)),3))];
    poss_early_fail_wins = [poss_early_fail_wins, i];
end
fail_temp = (poss_early_fail_win_vals(poss_early_fail_win_vals>min(poss_early_fail_win_vals)));
% fail_vals = min(fail_temp(fail_temp>min(fail_temp)));
fail_vals = min(fail_temp);
resp_early_fail_window = poss_early_fail_wins(find(poss_early_fail_win_vals == fail_vals)); %selects window with the peak response 
resp_early_fail_window = [resp_early_fail_window-1, resp_early_fail_window, resp_early_fail_window+1]; 

poss_late_fail_wins = [];  %using a sliding window to determine peak response to lever release 
poss_late_fail_win_vals = []; 
for i = pre_release_frames+1:pre_release_frames+1+round(300./double(ifi))
    poss_late_fail_win_vals = [poss_late_fail_win_vals; mean(mean(mean(fail_movie(:,:,i-1:i+1)),3))];
    poss_late_fail_wins = [poss_late_fail_wins, i];
end
fail_temp = (poss_late_fail_win_vals(poss_late_fail_win_vals>min(poss_late_fail_win_vals)));
% fail_vals = min(fail_temp(fail_temp>min(fail_temp)));
fail_vals = min(fail_temp);
resp_late_fail_window = poss_late_fail_wins(find(poss_late_fail_win_vals == fail_vals)); %selects window with the peak response 
resp_late_fail_window = [resp_late_fail_window-1, resp_late_fail_window, resp_late_fail_window+1]; 

poss_rel_wins = [];  %using a sliding window to determine peak response to lever release 
poss_rel_win_vals = []; 
for i = pre_release_frames+1:pre_release_frames+1+round(300./double(ifi))
    poss_rel_win_vals = [poss_rel_win_vals; mean(mean(mean(release_movie(:,:,i-1:i+1)),3))];
    poss_rel_wins = [poss_rel_wins, i];
end
resp_release_window = poss_rel_wins(find(poss_rel_win_vals == max(poss_rel_win_vals))); %selects window with the peak response 
resp_release_window = [resp_release_window-1, resp_release_window, resp_release_window+1]; 

poss_tf_wins = [];  %using a sliding window to determine peak response to lever release 
poss_tf_win_vals = []; 
for i = pre_release_frames+1:pre_release_frames+1+round(300./double(ifi))
    poss_tf_win_vals = [poss_tf_win_vals; mean(mean(mean(tooFastCorrect_movie(:,:,i-1:i+1)),3))];
    poss_tf_wins = [poss_tf_wins, i];
end
resp_tf_window = poss_tf_wins(find(poss_tf_win_vals == min(poss_tf_win_vals))); %selects window with the peak response 
resp_tf_window = [resp_tf_window-1, resp_tf_window, resp_tf_window+1]; 

% release
release_base = squeeze(mean(release_movie(:,:,base_release_window),3));
release_resp = squeeze(mean(release_movie(:,:,resp_release_window),3));
release_resp_avg = mean((release_resp-release_base),1);
release_resp_sem = std((release_resp-release_base),[],1)./sqrt(size(release_resp,1));

%success
success_base = squeeze(mean(success_movie(:,:,base_release_window),3));
% success_base = squeeze(mean(success_movie(:,:,pre_release_frames-5:pre_release_frames-2),3));
success_resp = squeeze(mean(success_movie(:,:,resp_success_window),3));
success_resp_avg = mean((success_resp-success_base),1);
success_resp_sem = std((success_resp-success_base),[],1)./sqrt(size(success_resp,1));
success_resp2 = squeeze(mean(success_movie(:,:,pre_release_frames+6:pre_release_frames+11),3));
success_resp1 = squeeze(mean(success_movie(:,:,pre_release_frames+4:pre_release_frames+5),3));
late_success_base = squeeze(mean(late_success_movie(:,:,base_release_window),3));
late_success_resp = squeeze(mean(late_success_movie(:,:,pre_release_frames-1:pre_release_frames+4),3));

%tooFast
tooFast_base = squeeze(mean(tooFastCorrect_movie(:,:,base_release_window),3));
%success_resp = squeeze(mean(tooFastCorrect_movie(:,:,resp_release_window),3));
tooFast_resp = squeeze(mean(tooFastCorrect_movie(:,:,resp_tf_window),3));
tooFast_resp2 = squeeze(mean(tooFastCorrect_movie(:,:,pre_release_frames+2:pre_release_frames+6),3));
tooFast_resp_avg = mean((tooFast_resp-tooFast_base),1);
tooFast_resp_sem = std((tooFast_resp-tooFast_base),[],1)./sqrt(size(tooFast_resp,1));


%fail
fail_base = squeeze(mean(fail_movie(:,:,base_release_window),3));
fail_resp = squeeze(mean(fail_movie(:,:,resp_fail_window),3));
fail_resp_avg = mean((fail_resp-fail_base),1);
fail_resp_sem = std((fail_resp-fail_base),[],1)./sqrt(size(fail_resp,1));

early_fail_base = squeeze(mean(early_fail_movie(:,:,base_release_window),3));
early_fail_resp = squeeze(mean(early_fail_movie(:,:,resp_late_fail_window),3));
late_fail_base = squeeze(mean(late_fail_movie(:,:,base_release_window),3));
late_fail_resp = squeeze(mean(late_fail_movie(:,:,resp_late_fail_window),3));

early_fail_base = squeeze(mean(early_fail_movie(:,:,base_release_window),3));
early_fail_resp1 = squeeze(mean(early_fail_movie(:,:,pre_release_frames:pre_release_frames+5),3));
early_fail_resp2 = squeeze(mean(early_fail_movie(:,:,pre_release_frames+6:pre_release_frames+11),3));
late_fail_base = squeeze(mean(late_fail_movie(:,:,base_release_window),3));
late_fail_resp1 = squeeze(mean(late_fail_movie(:,:,pre_release_frames+3:pre_release_frames+5),3));
late_fail_resp2 = squeeze(mean(late_fail_movie(:,:,pre_release_frames+6:pre_release_frames+11),3));
fail_resp1 = squeeze(mean(fail_movie(:,:,pre_release_frames+3:pre_release_frames+5),3));
fail_resp2 = squeeze(mean(fail_movie(:,:,pre_release_frames+6:pre_release_frames+11),3));

fake_cell = [];
if ~isempty(avg_late_fail)
    for i = 1:size(avg_late_fail,1)
        if mean(avg_late_fail(i,base_release_window_new)) > mean(avg_late_fail(i,resp_fail_window))*0.65
%         if mean(diff(avg_late_fail(i,base_release_window_new))) < 0.5*mean(diff(avg_late_fail(i,base_release_window_new(1)-2:base_release_window_new(1))))
            fake_cell = [fake_cell i];
            
        end
    end
end

fake_cell2 = [];
if ~isempty(avg_early_fail)
    for i = 1:size(avg_early_fail,1)
        if mean(avg_early_fail(i,base_release_window_new)) > mean(avg_early_fail(i,resp_fail_window))*0.65
%         if mean(diff(avg_late_fail(i,base_release_window_new))) < 0.5*mean(diff(avg_late_fail(i,base_release_window_new(1)-2:base_release_window_new(1))))
            fake_cell2 = [fake_cell2 i];
            
        end
    end
end

%by press
base_press_window = 1:pre_press_frames-round(300./double(ifi));
resp_press_window = pre_press_frames-round(100./double(ifi)):pre_press_frames+round(100./double(ifi));
[maxx, indp] = max(mean(press_long_movie(:,:,resp_press_window),1),[],3);
press_base = zeros(size(press_long_movie,1),nCells);
press_resp = zeros(size(press_long_movie,1),nCells);
for ic = 1:nCells
    press_base(:,ic) = squeeze(mean(press_long_movie(:,ic,base_press_window),3));
    press_resp(:,ic) = squeeze(mean(press_long_movie(:,ic,resp_press_window(indp(ic))-round(50./double(ifi)):resp_press_window(indp(ic))+round(50./double(ifi))),3));
end
press_resp_avg = mean((press_resp-press_base),1);
press_resp_sem = std((press_resp-press_base),[],1)./sqrt(size(press_resp,1));

%by licks
% if ~isempty(lick_data)
%     lick_amp_window1 = (pre_press_frames-round(200./double(ifi))):(pre_press_frames + round(200./double(ifi)));
%     lick_amp_window2 = (pre_press_frames-round(200./double(ifi))):(pre_press_frames + round(200./double(ifi)));
%     
%     [lickb_resp, lickb_base, lickb_resp_avg, lickb_resp_sem] = getPeak(lick_movie, lick_amp_window1, pre_press_frames, ifi, []);
%     [licks_resp, licks_base, licks_resp_avg, licks_resp_sem] = getPeak(single_lick_movie, lick_amp_window2, pre_press_frames, ifi, []);
% else
%     lickb_resp =[]; lickb_base =[];
%     licks_resp = []; licks_base =[];
%     
% end
if ~isempty(lick_data)
    pre_lick_frames = round(2500/double(ifi));%lick_data.buffer_lick;
    post_lick_frames = round(2500/double(ifi));
    base_lick_window = 1:pre_lick_frames-floor(500/double(ifi));
    poss_licks_wins = [];  %using a sliding window to determine peak response to lever release
    poss_licks_win_vals = [];
    for i = pre_lick_frames-round(300./double(ifi)):pre_lick_frames+1+round(300./double(ifi))
        poss_licks_win_vals = [poss_licks_win_vals; mean(mean(mean(single_lick_movie(:,:,i-1:i+1)),3))];
        poss_licks_wins = [poss_licks_wins, i];
    end
    resp_licks_window = poss_licks_wins(find(poss_licks_win_vals == max(poss_licks_win_vals))); %selects window with the peak response
    resp_licks_window = [resp_licks_window-3 : resp_licks_window+3]; %extended window 180529LG
    licks_base = squeeze(mean(single_lick_movie(:,:,base_lick_window),3));
    licks_resp = squeeze(mean(single_lick_movie(:,:,resp_licks_window),3));
    licks_resp_avg = mean((licks_resp-licks_base),1);
    licks_resp_sem = std((licks_resp-licks_base),[],1)./sqrt(size(licks_resp,1));
    
    poss_lickb_wins = [];  %using a sliding window to determine peak response to lever release
    poss_lickb_win_vals = [];
    for i = pre_lick_frames-round(300./double(ifi)):pre_lick_frames+1+round(300./double(ifi))
        poss_lickb_win_vals = [poss_lickb_win_vals; mean(mean(mean(lick_movie(:,:,i-1:i+1)),3))];
        poss_lickb_wins = [poss_lickb_wins, i];
    end
    resp_lickb_window = poss_lickb_wins(find(poss_lickb_win_vals == max(poss_lickb_win_vals))); %selects window with the peak response
    resp_lickb_window = [resp_lickb_window-3 : resp_lickb_window+3]; %extended window 180529LG
    lickb_base = squeeze(mean(lick_movie(:,:,base_lick_window),3));
    lickb_resp = squeeze(mean(lick_movie(:,:,resp_lickb_window),3));
    lickb_resp_avg = mean((lickb_resp-lickb_base),1);
    lickb_resp_sem = std((lickb_resp-lickb_base),[],1)./sqrt(size(lickb_resp,1));
    
%     if ~isempty(omitReward_movie) || ~isempty(itiReward_movie)
        nlickF = 3;
        tlick = ((-pre_lick_frames:nlickF:pre_lick_frames + 1.5).*double(ifi));
        [success_lick_freq, sem_success_lick] = binLicks(success_lick, ifi);
        
        tlick2 = ((-pre_lick_frames:nlickF:pre_lick_frames + 1.5).*double(ifi));
        [successlong_lick_freq, sem_successlong_lick] = binLicks(success_lick_long, ifi);
        [faillong_lick_freq, sem_faillong_lick] = binLicks(fail_lick_long, ifi);
%     end
    
    [late_fail_lick_freq, sem_late_fail_lick] = binLicks(late_fail_lick, double(ifi));
    [early_fail_lick_freq, sem_early_fail_lick] = binLicks(early_fail_lick, double(ifi));
    
    success_lick_corr = corrLick(success_movie, success_lick);
    fail_lick_corr = corrLick(fail_movie, fail_lick);
%     release_lick_corr = corrLick(cat(1, success_movie, fail_movie), [success_lick; fail_lick]);
    
else
    lickb_resp =[]; lickb_base =[]; success_lick_corr = [];
    licks_resp = []; licks_base =[]; fail_lick_corr = [];
end

if ~isempty(omitReward_movie)
    poss_OR_wins = [];  %using a sliding window to determine peak response to lever release
    poss_OR_win_vals = [];
    for i = pre_release_frames+1:pre_release_frames+1+round(300./double(ifi))
        poss_OR_win_vals = [poss_OR_win_vals; mean(mean(mean(omitReward_movie(:,:,i-1:i+1)),3))];
        poss_OR_wins = [poss_OR_wins, i];
    end
    resp_OR_window = poss_OR_wins(find(poss_OR_win_vals == max(poss_OR_win_vals))); %selects window with the peak response
    resp_OR_window = [resp_OR_window-1, resp_OR_window, resp_OR_window+1];
    OR_base = squeeze(mean(omitReward_movie(:,:,base_release_window),3));
    OR_resp = squeeze(mean(omitReward_movie(:,:,resp_OR_window),3));
    OR_resp2 = squeeze(mean(omitReward_movie(:,:,pre_release_frames+6:pre_release_frames+11),3)); % same window as late earlys
    OR_resp1 = squeeze(mean(omitReward_movie(:,:,pre_release_frames+4:pre_release_frames+5),3)); % same window as early earlys
    OR_resp_avg = mean((OR_resp-OR_base),1);
    OR_resp_sem = std((OR_resp-OR_base),[],1)./sqrt(size(OR_resp,1));
    
%     succ_OR_resp = squeeze(mean(success_movie(:,:,pre_release_frames+6:pre_release_frames+11),3));
    
    [omitRlong_lick_freq, sem_omitRlong_lick] = binLicks(omitR_lick_long, ifi);
    
    [omitR_lick_freq, sem_omitR_lick] = binLicks(omitR_lick, ifi);
    
    omitR_lick_corr = corrLick(omitReward_movie, omitR_lick);
    
    
else
    OR_base = [];
    OR_resp = []; omitR_lick_corr=[]; OR_resp1 = []; OR_resp2 = [];
end

if ~isempty(itiReward_movie)
    poss_IR_wins = [];  %using a sliding window to determine peak response to lever release
    poss_IR_win_vals = [];
    for i = pre_release_frames+1:pre_release_frames+1+round(300./double(ifi))
        poss_IR_win_vals = [poss_IR_win_vals; mean(mean(mean(itiReward_movie(:,:,i-1:i+1)),3))];
        poss_IR_wins = [poss_IR_wins, i];
    end
    resp_IR_window = poss_IR_wins(find(poss_IR_win_vals == max(poss_IR_win_vals))); %selects window with the peak response
    resp_IR_window = [resp_IR_window-1, resp_IR_window, resp_IR_window+1];
    IR_base = squeeze(mean(itiReward_movie(:,:,base_release_window),3));
    IR_resp = squeeze(mean(itiReward_movie(:,:,resp_IR_window),3));
    IR_resp_avg = mean((IR_resp-IR_base),1);
    IR_resp_sem = std((IR_resp-IR_base),[],1)./sqrt(size(IR_resp,1));
    
    
    [itiR_lick_freq, sem_itiR_lick] = binLicks(itiR_lick, ifi);
    [itiRlong_lick_freq, sem_itiRlong_lick] = binLicks(itiR_lick_long, ifi);
    
    itiR_lick_corr = corrLick(itiReward_movie, itiR_lick);
else
    IR_base = [];
    IR_resp = [];
    itiR_lick_corr =[];
end

%3.b find cell response amplitude
% press_amp_window = (pre_press_frames-round(200./double(ifi))):(pre_press_frames + round(200./double(ifi)));
% default_p  = round(200./double(ifi)) + 2;
% [press_resp, press_base, press_resp_avg, press_resp_sem] = getPeak(press_long_movie, press_amp_window, pre_press_frames, ifi, default_p);
% 
% release_amp_window = (pre_release_frames-round(200./double(ifi))):(pre_release_frames + round(200./double(ifi)));
% default_p  = 3;
% [release_resp, release_base, release_resp_avg, release_resp_sem] = getPeak(release_long_movie, release_amp_window, pre_release_frames, ifi, default_p);
% [success_resp, success_base, success_resp_avg, success_resp_sem] = getPeak(success_movie, release_amp_window, pre_release_frames, ifi, default_p);
% [tooFast_resp, tooFast_base, tooFast_resp_avg, tooFast_resp_sem] = getPeak(tooFastCorrect_movie, release_amp_window, pre_release_frames, ifi, default_p);
% [fail_resp, fail_base, fail_resp_avg, fail_resp_sem] = getPeak(fail_movie, release_amp_window, pre_release_frames, ifi, default_p);


save([dest_sub '_cell_resp_LG.mat'], 'press_base', 'press_resp', 'success_base', 'success_resp', 'tooFast_resp', 'tooFast_base', 'fail_base', ...
    'fail_resp', 'early_fail_resp','late_fail_resp', 'early_fail_base', 'late_fail_base', 'release_resp', 'release_base', 'lickb_resp', 'lickb_base', 'licks_resp', 'licks_base', 'OR_base', 'OR_resp','IR_base', 'IR_resp',...
    'fail_resp1', 'fail_resp2', 'success_resp1', 'late_fail_resp1', 'late_fail_resp2', 'late_fail_base', 'fake_cell', 'early_fail_resp1', 'early_fail_resp2','early_fail_base', 'tooFast_resp2',...
    'late_success_base', 'late_success_resp', 'OR_resp1', 'OR_resp2', 'success_resp2');

save([dest_sub '_dff_lick_LG.mat'], 'success_lick_corr', 'fail_lick_corr', 'itiR_lick_corr', 'omitR_lick_corr');

%% 2. ttest for significant responses
[release_h, release_p] = ttest(release_base, release_resp, 'dim', 1, 'tail', 'both');
[success_h, success_p] = ttest(success_base, success_resp, 'dim', 1, 'tail', 'both');
[tooFast_h, tooFast_p] = ttest(tooFast_base, tooFast_resp, 'dim', 1, 'tail', 'both');
[fail_h, fail_p] = ttest(fail_base, fail_resp, 'dim', 1, 'tail', 'both');
[press_h, press_p] = ttest(press_base, press_resp, 'dim', 1, 'tail', 'both');
if ~isempty(lick_data)
    [lickb_h, lickb_p] = ttest(lickb_base, lickb_resp, 'dim', 1, 'tail', 'both');
    [licks_h, licks_p] = ttest(licks_base, licks_resp, 'dim', 1, 'tail', 'both');
end

[release_h, release_p] = ttest(release_base, release_resp, 'dim', 1, 'tail', 'left');
[success_h, success_p] = ttest(success_base, success_resp, 'dim', 1, 'tail', 'left');
[tooFast_h, tooFast_p] = ttest(tooFast_base, tooFast_resp, 'dim', 1, 'tail', 'left');
[fail_h, fail_p] = ttest(fail_base, fail_resp, 'dim', 1, 'tail', 'left');
[press_h, press_p] = ttest(press_base, press_resp, 'dim', 1, 'tail', 'left');
if ~isempty(lickb_resp)
    [lickb_h, lickb_p] = ttest(lickb_base, lickb_resp, 'dim', 1, 'tail', 'left');
else
    lickb_h = nan(1,size(lickb_resp,2));
end
if ~isempty(licks_resp)
    [licks_h, licks_p] = ttest(licks_base, licks_resp, 'dim', 1, 'tail', 'left');
    [licks_eh, licks_ep] = ttest(licks_base, licks_resp, 'dim', 2, 'tail', 'left');
    licks_resp_per = sum(licks_eh)/length(licks_eh);
else
    licks_h = nan(1,size(licks_resp,2));
    licks_eh = nan(1,size(licks_resp,2));
end

if ~isempty(OR_base)
    [OR_h, OR_p] = ttest(OR_base, OR_resp, 'dim', 1, 'tail', 'left');
else
    OR_h = [];
end

if ~isempty(IR_base)
    [IR_h, IR_p] = ttest(IR_base, IR_resp, 'dim', 1, 'tail', 'left');
else
    IR_h = [];
end


%% 3. define cells by response to event
release_resp_cells = find(release_h);
success_resp_cells = find(success_h);
tooFast_resp_cells = find(tooFast_h);
fail_resp_cells = find(fail_h);
press_resp_cells = find(press_h);
if ~isempty(lick_data)
    lickb_resp_cells = find(lickb_h);
    licks_resp_cells = find(licks_h);
else
    lickb_resp_cells =[]; licks_resp_cells = [];
end

if ~isempty(OR_h)
    OR_resp_cells = find(OR_h);
else
    OR_resp_cells = [];
end
if ~isempty(IR_h)
    IR_resp_cells = find(IR_h);
else
    IR_resp_cells = [];
end

fail_only_cells = fail_resp_cells(find(ismember(fail_resp_cells, success_resp_cells)==0));
success_only_cells = success_resp_cells(find(ismember(success_resp_cells, fail_resp_cells)==0));
fail_and_success_cells = find(and(fail_h, success_h));
press_and_success_cells = find(and(success_h,press_h));
press_notsuccess_cells = press_resp_cells(find(ismember(press_resp_cells, success_resp_cells)==0));
success_notpress_cells = success_resp_cells(find(ismember(success_resp_cells, press_resp_cells)==0));
noresponse_cells = find((success_h+fail_h+press_h)==0);

if ~isempty(lick_data)
    %lickbs_h = or(lickb_h, licks_h);
    lick_cells = unique([find(lickb_h==1) find(licks_h==1)]);
    notlick_cells = intersect(find(lickb_h==0), find(licks_h==0));
%     notlick_rel_cells = intersect(notlick_cells, release_resp_cells);
%     lick_notsuccess_cells = lick_cells(find(ismember(lick_cells, success_resp_cells)==0));
%     lick_success_cells = find(and(lickbs_h, success_h));
%     licks_success_cells = find(and(licks_h, success_h));
%     lickb_success_cells = find(and(lickb_h, success_h));
%     success_notlick_cells = success_resp_cells(find(ismember(success_resp_cells, lick_cells)==0));
else
    lick_cells = []; notlick_cells = []; notlick_rel_cells = []; lick_success_cells = [];
    licks_success_cells = []; lickb_success_cells = []; success_notlick_cells = []; lick_notsuccess_cells = [];
end

% find transient and sustain cells
success_TC_all_RL = avg_success(release_resp_cells,1:(pre_release_frames*2+1));
[~,indx] = max(success_TC_all_RL,[],2);
tot_cell = size(success_TC_all_RL,1);
frames = size(success_TC_all_RL,2);
trans_cell = zeros(tot_cell,1);
for ic = 1:tot_cell
   
    if frames - indx(ic) < 2
        
        temp_tc = success_TC_all_RL(ic,:);
        [~,indx(ic)] = max(diff(temp_tc));
        indx(ic) = indx(ic) + 2;
    end
    if frames - indx(ic) > 2
        temp_tc = success_TC_all_RL(ic,indx(ic):end);
        if mean(temp_tc)/max(temp_tc) < 0.5
            trans_cell(ic) = 1;
        else
            trans_cell(ic) = 0.5;
        end
            
    end
end
trans_cell_ind = find(trans_cell==1); sustain_cell_ind = find(trans_cell==0.5);

save([dest_sub '_cell_categories_LG.mat'], 'release_resp_cells', 'success_resp_cells', 'tooFast_resp_cells', 'fail_resp_cells', 'press_resp_cells',...
    'noresponse_cells', 'fail_only_cells', 'success_only_cells', 'fail_and_success_cells', 'press_and_success_cells', 'press_notsuccess_cells', ...
    'success_notpress_cells','lickb_resp_cells', 'licks_resp_cells', 'trans_cell_ind', 'sustain_cell_ind', 'lick_cells', 'notlick_cells', 'OR_resp_cells', 'IR_resp_cells');

% % %test success vs fail responses (in cells that respond to both)
% % [fail_v_success_h, fail_v_success_p] = ttest(fail_resp_avg(1,fail_and_success_cells),success_resp_avg(1,fail_and_success_cells), 'tail', 'left');
% % save([dest_sub '_pvals.mat'], 'fail_v_success_p', 'release_p','success_p', 'fail_p', 'press_p', 'tooFast_p', 'release_h', 'success_h', 'fail_h', 'press_h', 'tooFast_h');
% % 
% % 
% % 
% % 
% % tt =((-pre_release_frames:post_release_frames).*double(ifi));
% % % %overlay average success/failure for each ROI
% % fig=figure;
% % avg_all = [avg_early_fail avg_late_fail];
% % ymax = max(max(avg_all,[],2),[],1);
% % ymin = min(min(avg_all,[],2),[],1);
% % n_h = avg_all;
% % for ic = 1:size(n_h,1)
% %     subplot(n,n2,ic)
% %     hold on
% %     plot(tt,avg_early_fail(ic,:),'k')
% %     plot(tt,avg_late_fail(ic,:),'r')
% %     ylim([ymin*1.1 ymax*1.1])
% %     xlim([tt(1) tt(end)])
% % end
% % supertitle([date ' ' mouse ' Average release: Early- m (n = ' num2str(size(success_movie,1)) ' trials); Late- black (n = ' num2str(size(fail_movie,1)) ' trials)'])
% % orient landscape
% % % saveas(fig, [dest_sub '_earlylate_fail_allTCs.fig']);
% % % % print([dest_sub '_earlylate_fail_allTCs.eps'], '-depsc');
% % % % print([dest_sub '_earlylate_fail_allTCs.pdf'], '-dpdf');
% % release_h1 = release_h;
% % release_h1(fake_cell2) = 0;
% % 
% % release_h = 0*release_h;
% % release_h(fake_cell2(21:39)) = 0;
% % tt =((-pre_release_frames:post_release_frames).*double(ifi));
% % fig = figure;
% % hold on;
% % avg_early_fail_temp = avg_early_fail;
% % avg_early_fail_temp(:,1:end-2) = avg_early_fail(:,3:end);
% % avg_early_fail_temp(:,end-2:end) = avg_early_fail(:,1:3);
% % shadedErrorBar(tt,mean(avg_early_fail_temp(find(release_h1),:),1)*0.4,std(avg_early_fail_temp(find(release_h1),:),[],1)./sqrt(sum(release_h,2)),'r')
% % shadedErrorBar(tt,mean(avg_late_fail(find(release_h),:),1),std(avg_late_fail(find(release_h),:),[],1)*0.9./sqrt(sum(release_h,2)),'k')
% % 
% % fig = figure;
% % hold on;
% % shadedErrorBar(tt,mean(avg_early_fail,1)+0.02,std(avg_early_fail,[],1)./sqrt(size(avg_early_fail,1)),'r')
% % shadedErrorBar(tt,mean(avg_late_fail,1),std(avg_late_fail,[],1)*0.8./sqrt(size(avg_late_fail,1)),'k')
% % 
% % title(['Avg success resp p<0.05; n = ' num2str(sum(release_h,2))])
% % ylabel('dF/F')
% % xlabel('Time (ms)')
% % supertitle([date ' ' mouse ' Average release: Late- black; Early- red'])
% % saveas(fig, [dest_sub '_earlylate_fail_avgTCs.fig']);
% % print([dest_sub '_earlylate_fail_avgTCs.eps'], '-depsc');
% % print([dest_sub '_earlylate_fail_avgTCs.pdf'], '-dpdf');
% 
% 
% % fig = figure;
% % tt =((-pre_release_frames:post_release_frames).*double(ifi));
% % shadedErrorBar(tt,mean(avg_early_fail(find(release_h),:),1),std(avg_early_fail(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'m')
% % hold on
% % shadedErrorBar(tt,mean(avg_late_fail(find(release_h),:),1),std(avg_late_fail(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'r')
% % bar(tlick(1:end-1), (late_fail_lick_freq/500), 'r');
% % bar(tlick(1:end-1), (early_fail_lick_freq)/500, 'm'); 
% % errorbar(tlick(1:end-1), (early_fail_lick_freq)/500, sem_early_fail_lick/500, '.b');
% % errorbar(tlick(1:end-1), (late_fail_lick_freq)/500, sem_late_fail_lick/500, '.r');
% % % xlim([-500 500])
% % ylabel('dF/F')
% % xlabel('Time (ms)')
% % title(['Avg resp p<0.05; n = ' num2str(sum(release_h,2)) ' early: magenta, late early: red'])
% % title(['early: magenta, late early: red'])
% % saveas(fig, [dest_sub '_earlylate_fail_lick_avgTCs.fig']);
% % print([dest_sub '_earlylate_fail_lick_avgTCs.eps'], '-depsc');
% % print([dest_sub '_earlylate_fail_lick_avgTCs.pdf'], '-dpdf');
% 
% 
% % %% 4. plotting
% % tt =((-pre_release_frames:post_release_frames).*double(ifi));
% % %overlay average success/failure for each ROI
% % fig=figure;
% % avg_all = [avg_success avg_fail];
% % ymax = max(max(avg_all,[],2),[],1);
% % ymin = min(min(avg_all,[],2),[],1);
% % for ic = 1:nCells
% %     subplot(n,n2,ic)
% %     errorbar(tt,avg_success(ic,:), sem_success(ic,:),'k')
% %     hold on;
% %     errorbar(tt,avg_fail(ic,:), sem_fail(ic,:),'r')
% %     ylim([ymin*1.1 ymax*1.1])
% %     xlim([tt(1) tt(end)])
% % end
% % supertitle([date ' ' mouse ' Average release: Success- black (n = ' num2str(size(success_movie,1)) ' trials); Failure- red (n = ' num2str(size(fail_movie,1)) ' trials)'])
% % orient landscape
% % saveas(fig, [dest_sub '_success_fail_allTCs.fig']);
% % print([dest_sub '_success_fail_allTCs.eps'], '-depsc');
% % print([dest_sub '_success_fail_allTCs.pdf'], '-dpdf');
% % 
% % % overlay average success/failure across ROIs
% % fig=figure; 
% % subplot(2,2,1)
% % % errorbar(tt,avg_success_all, sem_success_all,'k')
% % hold on;
% % % errorbar(tt,avg_fail_all, sem_fail_all,'r')
% % errorbar(tt,avg_lick_all, sem_lick_all,'b')
% % errorbar(tt,avg_success_all, sem_success_all,'k')
% % errorbar(tt,avg_single_lick_all, sem_single_lick_all,'g')
% % errorbar(tt,avg_fail_all, sem_fail_all,'r')
% % title(['Avg all cells; n = ' num2str(size(avg_success,1))])
% % ylabel('dF/F')
% % xlabel('Time (ms)')
% % subplot(2,2,2)
% % % errorbar(tt,mean(avg_success(find(release_h),:),1),std(avg_success(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'k')
% % hold on;
% % % errorbar(tt,mean(avg_fail(find(release_h),:),1),std(avg_fail(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'r')
% % errorbar(tt,mean(avg_lick(find(lickb_h),:),1),std(avg_lick(find(lickb_h),:),[],1)./sqrt(sum(lickb_h,2)),'b')
% % errorbar(tt,mean(avg_single_lick(find(licks_h),:),1),std(avg_single_lick(find(licks_h),:),[],1)./sqrt(sum(licks_h,2)),'g')
% % errorbar(tt,mean(avg_success(find(success_h),:),1),std(avg_success(find(success_h),:),[],1)./sqrt(sum(success_h,2)),'k')
% % errorbar(tt,mean(avg_fail(find(fail_h),:),1),std(avg_fail(find(fail_h),:),[],1)./sqrt(sum(fail_h,2)),'r')
% % 
% % title(['Avg resp p<0.05; n = ' num2str(sum(release_h,2))])
% % ylabel('dF/F')
% % xlabel('Time (ms)')
% % subplot(2,2,3)
% % % errorbar(tt,mean(avg_success(find(success_h),:),1),std(avg_success(find(success_h),:),[],1)./sqrt(sum(success_h,2)),'k')
% % hold on;
% % % errorbar(tt,mean(avg_fail(find(success_h),:),1),std(avg_fail(find(success_h),:),[],1)./sqrt(sum(success_h,2)),'r')
% % errorbar(tt,mean(avg_lick(find(success_h),:),1),std(avg_lick(find(success_h),:),[],1)./sqrt(sum(success_h,2)),'b')
% % errorbar(tt,mean(avg_single_lick(find(success_h),:),1),std(avg_single_lick(find(success_h),:),[],1)./sqrt(sum(success_h,2)),'g')
% % errorbar(tt,mean(avg_success(find(success_h),:),1),std(avg_success(find(success_h),:),[],1)./sqrt(sum(success_h,2)),'k')
% % errorbar(tt,mean(avg_fail(find(success_h),:),1),std(avg_fail(find(success_h),:),[],1)./sqrt(sum(success_h,2)),'r')
% % 
% % title(['Avg success resp p<0.05; n = ' num2str(sum(success_h,2))])
% % ylabel('dF/F')
% % xlabel('Time (ms)')
% % supertitle([date ' ' mouse ' Average release: Success- black; Failure- red Lickbout- blue Lick- green'])
% % saveas(fig, [dest_sub '_success_fail_avgTCs.fig']);
% % print([dest_sub '_success_fail_avgTCs.eps'], '-depsc');
% % print([dest_sub '_success_fail_avgTCs.pdf'], '-dpdf');
% % 
% % % overlay average omit Reward and normal reward
% % 
% % if ~isempty(OR_base)
% %     fig=figure;
% %     tt2 = ((-pre_release_frames2:post_release_frames2).*double(ifi));
% %     avg_all = [avg_success_long avg_OR_long];
% %     ymax = max(max(avg_all,[],2),[],1);
% %     ymin = min(min(avg_all,[],2),[],1);
% %     for ic = 1:nCells
% %         subplot(n,n2,ic)
% %         hold on;
% %         errorbar(tt2,avg_OR_long(ic,:), sem_OR_long(ic,:),'r')
% %         errorbar(tt2,avg_success_long(ic,:), sem_success_long(ic,:),'k')
% %         
% %         ylim([ymin*1.1 ymax*1.1])
% %         xlim([tt2(1) tt2(end)])
% %     end
% %     supertitle([date ' ' mouse ' Average release: Success- black (n = ' num2str(size(success_movie,1)) ' trials); Omit Reward- red (n = ' num2str(size(omitReward_movie,1)) ' trials)'])
% %     orient landscape
% % %     saveas(fig, [dest_sub '_success_omitR_allTCs.fig']);
% % %     print([dest_sub '_success_omitR_allTCs.eps'], '-depsc');
% % %     print([dest_sub '_success_omitR_allTCs.pdf'], '-dpdf');
% % %     
% %     fig = figure;
% %    
% %     shadedErrorBar(tt2,mean(avg_success_long(find(release_h),:),1)+0.023,std(avg_success_long(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'k')
% %     hold on
% %     shadedErrorBar(tt2,mean(avg_OR_long(find(release_h),:),1)+0.01,std(avg_OR_long(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'r')
% %     bar(tlick2(1:end-1), successlong_lick_freq/500, 'k'); errorbar(tlick2(1:end-1), successlong_lick_freq/500, sem_successlong_lick/500, '.k');
% %     bar(tlick2(1:end-1), omitRlong_lick_freq/500, 'r'); errorbar(tlick2(1:end-1), omitRlong_lick_freq/500, sem_omitRlong_lick/500, '.r');
% %     
% %     ylabel('dF/F')
% %     xlabel('Time (ms)')
% %     title(['Avg resp p<0.05; n = ' num2str(sum(release_h,2)) ' Success: black, Omit Reward: red'])
% %     saveas(fig, [dest_sub '_success_omitR_avgTCs.fig']);
% %     print([dest_sub '_success_omitR_avgTCs.eps'], '-depsc');
% %     print([dest_sub '_success_omitR_avgTCs.pdf'], '-dpdf');
% %     
% %     
% %     fig = figure;
% %   
% %     shadedErrorBar(tt2,mean(avg_fail_long(find(release_h),:),1),std(avg_fail_long(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'b')
% %     hold on
% %     shadedErrorBar(tt2,mean(avg_OR_long(find(release_h),:),1),std(avg_OR_long(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'r')
% %     bar(tlick2(1:end-1), omitRlong_lick_freq/500, 'r'); errorbar(tlick2(1:end-1), omitRlong_lick_freq/500, sem_omitRlong_lick/500, '.r');
% %     bar(tlick2(1:end-1), faillong_lick_freq/500, 'b'); errorbar(tlick2(1:end-1), faillong_lick_freq/500, sem_faillong_lick/500, '.b');
% %     
% %     
% %     ylabel('dF/F')
% %     xlabel('Time (ms)')
% %     title(['Avg resp p<0.05; n = ' num2str(sum(release_h,2)) ' Early: blue, Omit Reward: red'])
% %     saveas(fig, [dest_sub '_fail_omitR_avgTCs.fig']);
% %     print([dest_sub '_fail_omitR_avgTCs.eps'], '-depsc');
% %     print([dest_sub '_fail_omitR_avgTCs.pdf'], '-dpdf');
% %     
% %     omitRewardIndx = find(trial_outcome.omitRewardIndx);
% %     prev_trial = []; omitR_prevR_TC = []; omitR_prevR_lick = [];
% %     omitR_prevNR_TC = []; omitR_prevNR_lick = [];
% %     
% %     for j = 1:size(omitReward_movie,1)
% %         if trial_outcome.hasReward(omitRewardIndx(j) - 1) == 1
% %             omitR_prevR_TC = squeeze(omitReward_movie(j, :, :));
% %             omitR_prevR_lick = [omitR_prevR_lick; omitR_lick(j,:)];
% %         else
% %             omitR_prevNR_TC = squeeze(omitReward_movie(j, :, :));
% %             omitR_prevNR_lick = [omitR_prevNR_lick; omitR_lick(j,:)];
% %         end
% %     end
% %     
% %     
% %     fig = figure;
% %     hold on
% %     errorbar(tt,mean(omitR_prevR_TC(find(release_h),:),1),std(omitR_prevR_TC(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'b')
% %     errorbar(tt,mean(omitR_prevNR_TC(find(release_h),:),1),std(omitR_prevNR_TC(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'r')
% %     ylabel('dF/F')
% %     xlabel('Time (ms)')
% %     title(['Avg resp p<0.05; n = ' num2str(sum(release_h,2)) ' Omit after reward: blue, Omit after no reward: red'])
% %     saveas(fig, [dest_sub '_omitR_previousTrial_avgTCs.fig']);
% %     print([dest_sub '_omitR_previousTrial_avgTCs.eps'], '-depsc');
% %     print([dest_sub '_omitR_previousTrial_avgTCs.pdf'], '-dpdf');
% %     
% %     save([dest_sub '_omitR_previousTrial.mat'], 'omitR_prevR_TC', 'omitR_prevNR_TC', 'omitR_prevR_lick', 'omitR_prevNR_lick');
% % end
% % 
% % 
% % if ~isempty(IR_base)
% %     fig=figure;
% %     
% %     avg_all = [avg_success avg_IR];
% %     ymax = max(max(avg_all,[],2),[],1);
% %     ymin = min(min(avg_all,[],2),[],1);
% %     for ic = 1:nCells
% %         subplot(n,n2,ic)
% %         errorbar(tt,avg_success(ic,:), sem_success(ic,:),'k')
% %         hold on;
% %         errorbar(tt,avg_IR(ic,:), sem_IR(ic,:),'g')
% %         ylim([ymin*1.1 ymax*1.1])
% %         xlim([tt(1) tt(end)])
% %     end
% %     supertitle([date ' ' mouse ' Average release: Success- black (n = ' num2str(size(success_movie,1)) ' trials); Unexp Reward- green (n = ' num2str(size(itiReward_movie,1)) ' trials)'])
% %     orient landscape
% %     saveas(fig, [dest_sub '_success_itiR_allTCs.fig']);
% %     print([dest_sub '_success_itiR_allTCs.eps'], '-depsc');
% %     print([dest_sub '_success_itiR_allTCs.pdf'], '-dpdf');
% %     
% %     tt2 = ((-pre_release_frames2:post_release_frames2).*double(ifi));
% %     fig = figure;
% %     errorbar(tt2,mean(avg_success_long(find(release_h),:),1),std(avg_success_long(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'k')
% %     hold on
% %     errorbar(tt2,mean(avg_IR_long(find(release_h),:),1),std(avg_IR_long(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'g')
% %     bar(tlick2(1:end-1), successlong_lick_freq/500, 'k'); errorbar(tlick2(1:end-1), successlong_lick_freq/500, sem_successlong_lick/500, '.k');
% %     bar(tlick2(1:end-1), itiRlong_lick_freq/500, 'g'); errorbar(tlick2(1:end-1), itiRlong_lick_freq/500, sem_itiRlong_lick/500, '.g');
% %     
% %     ylabel('dF/F')
% %     xlabel('Time (ms)')
% %     title(['Avg resp p<0.05; n = ' num2str(sum(release_h,2)) ' Success: black, iti Unexpected Reward: green'])
% %     saveas(fig, [dest_sub '_success_itiR_avgTCs.fig']);
% %     print([dest_sub '_success_itiR_avgTCs.eps'], '-depsc');
% %     print([dest_sub '_success_itiR_avgTCs.pdf'], '-dpdf');
% % end
% % 
% % %overlay average press/release for each ROI
% % fig=figure;
% % avg_all = [avg_press avg_release];
% % ymax = max(max(avg_all,[],2),[],1);
% % ymin = min(min(avg_all,[],2),[],1);
% % for ic = 1:nCells
% %     subplot(n,n2,ic)
% %     errorbar(tt,avg_press(ic,:), sem_press(ic,:),'c')
% %     hold on;
% %     errorbar(tt,avg_release(ic,:), sem_release(ic,:),'k')
% %     ylim([ymin*1.1 ymax*1.1])
% %     xlim([tt(1) tt(end)])
% % end
% % supertitle([date ' ' mouse ' Average release: black (n = ' num2str(size(release_movie,1)) ' trials); Press- cyan (n = ' num2str(size(press_long_movie,1)) ' trials)'])
% % orient landscape
% % saveas(fig, [dest_sub '_press_release_allTCs.fig']);
% % print([dest_sub '_press_release_allTCs.eps'], '-depsc');
% % print([dest_sub '_press_release_allTCs.pdf'], '-dpdf');
% % %overlay average press/release across ROIs
% % fig=figure;
% % subplot(2,2,1)
% % errorbar(tt,avg_release_all, sem_release_all,'k')
% % hold on;
% % errorbar(tt,avg_press_all, sem_press_all,'c')
% % title(['Avg all cells; n = ' num2str(size(avg_release,1))])
% % ylabel('dF/F')
% % xlabel('Time (ms)')
% % subplot(2,2,2)
% % errorbar(tt+1500,mean(avg_release(find(release_h),:),1),std(avg_release(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'k')
% % hold on;
% % errorbar(tt,mean(avg_press(find(release_h),:),1),std(avg_press(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'c')
% % errorbar(tt+1500,mean(avg_fail(find(release_h),:),1),std(avg_fail(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'r')
% % vline([0 1], '--k'); vline(1500, '--k');
% % set(gca, 'XTick', [0, 1500], 'XTickLabel', {'press', 'release, early'});
% % title(['Avg release resp p<0.05; n = ' num2str(sum(release_h,2))])
% % ylabel('dF/F')
% % xlabel('Time (ms)')
% % subplot(2,2,3)
% % errorbar(tt,mean(avg_release(find(press_h),:),1),std(avg_release(find(press_h),:),[],1)./sqrt(sum(press_h,2)),'k')
% % hold on;
% % errorbar(tt,mean(avg_press(find(press_h),:),1),std(avg_press(find(press_h),:),[],1)./sqrt(sum(press_h,2)),'c')
% % title(['Avg press resp p<0.05; n = ' num2str(sum(press_h,2))])
% % ylabel('dF/F')
% % xlabel('Time (ms)')
% % supertitle([date ' ' mouse ' Average release: Release- black; Press- cyan'])
% % saveas(fig, [dest_sub '_press_release_avgTCs.fig']);
% % print([dest_sub '_press_release_avgTCs.eps'], '-depsc');
% % print([dest_sub '_press_release_avgTCs.pdf'], '-dpdf');
% % 
% % %overlay all lick bout/lick for each ROI
% % if ~isempty(lick_data)
% %     fig=figure;
% %     avg_all = [avg_lick avg_single_lick];
% %     ymax = max(max(avg_all,[],2),[],1);
% %     ymin = min(min(avg_all,[],2),[],1);
% %     for ic = 1:nCells
% %         subplot(n,n2,ic)
% %         errorbar(tt,avg_lick(ic,:), sem_lick(ic,:),'b')
% %         hold on;
% %         errorbar(tt,avg_single_lick(ic,:), sem_single_lick(ic,:),'g')
% %         ylim([ymin*1.1 ymax*1.1])
% %         xlim([tt(1) tt(end)])
% %     end
% %     supertitle([date ' ' mouse ' Average lick bout: blue (n = ' num2str(size(lick_movie,1)) ' trials); single lick- green (n = ' num2str(size(single_lick_movie,...
% %         1)) ' trials)'])
% %     orient landscape
% %     saveas(fig, [dest_sub '_lickbout_singlelick_allTCs.fig']);
% %     print([dest_sub '_lickbout_singlelick_allTCs.eps'], '-depsc');
% %     print([dest_sub '_lickbout_singlelick_allTCs.pdf'], '-dpdf');
% %     
% %     fig=figure;
% %     avg_all = [avg_single_lick avg_lick avg_success];
% %     ymax = max(max(avg_all,[],2),[],1);
% %     ymin = min(min(avg_all,[],2),[],1);
% %     
% %     if length(lick_notsuccess_cells) <= 4
% %         ns = 2;
% %     else
% %         ns = ceil(sqrt(length(lick_notsuccess_cells)));
% %     end
% %     for ic = 1:length(lick_notsuccess_cells)
% %         n_ic = lick_notsuccess_cells(ic);
% %         subplot(ns,ceil(length(lick_notsuccess_cells)/ns),ic)
% %         if ~isempty(find(licks_resp_cells==n_ic))
% %             errorbar(tt,avg_single_lick(n_ic,:), sem_single_lick(n_ic,:),'b')
% %         else
% %             errorbar(tt,avg_lick(n_ic,:), sem_single_lick(n_ic,:),'g')
% %         end
% %         hold on;
% %         errorbar(tt,avg_success(n_ic,:), sem_success(n_ic,:),'k')
% %         ylim([ymin*1.1 ymax*1.1])
% %         xlim([tt(1) tt(end)])
% %     end
% %     supertitle([date ' ' mouse ' Average single lick: blue (n = ' num2str(size(single_lick_movie,1)) ' trials);' ' lickbout: green (n= ' num2str(size(lick_movie,1)) ' lick not success resp- black (n = ' num2str(size(success_movie,...
% %         1)) ' trials) nCells = ' num2str(length(lick_notsuccess_cells))])
% %     orient landscape
% %     saveas(fig, [dest_sub '_lickNotSuccess_allTCs.fig']);
% %     print([dest_sub '_lickNotSuccess_allTCs.eps'], '-depsc');
% %     print([dest_sub '_lickNotSuccess_allTCs.pdf'], '-dpdf');
% %     
% %     fig=figure;
% %     avg_all = [avg_single_lick avg_lick avg_success];
% %     ymax = max(max(avg_all,[],2),[],1);
% %     ymin = min(min(avg_all,[],2),[],1);
% %     if length(lick_success_cells) <= 4
% %         ns = 2;
% %     else
% %         ns = ceil(sqrt(length(lick_success_cells)));
% %     end
% %     for ic = 1:length(lick_success_cells)
% %         n_ic = lick_success_cells(ic);
% %         subplot(ns,ceil(length(lick_success_cells)/ns),ic)
% %         if ~isempty(find(licks_resp_cells==n_ic))
% %             errorbar(tt,avg_single_lick(n_ic,:), sem_single_lick(n_ic,:),'b')
% %         else
% %             errorbar(tt,avg_lick(n_ic,:), sem_single_lick(n_ic,:),'g')
% %         end
% %         hold on;
% %         errorbar(tt,avg_success(n_ic,:), sem_success(n_ic,:),'k')
% %         ylim([ymin*1.1 ymax*1.1])
% %         xlim([tt(1) tt(end)])
% %     end
% %     supertitle([date ' ' mouse ' Average single lick: blue (n = ' num2str(size(single_lick_movie,1)) 'trials);' ' lickbout: green (n= ' num2str(size(lick_movie,1)) 'lick not success resp- black (n = ' ,...
% %         num2str(length(lick_success_cells))])
% %     orient landscape
% %     saveas(fig, [dest_sub '_lickSuccess_allTCs.fig']);
% %     print([dest_sub '_lickSuccess_allTCs.eps'], '-depsc');
% %     print([dest_sub '_lickSuccess_allTCs.pdf'], '-dpdf');
% %     
% %     fig = figure;
% %     hold on
% %     errorbar(tt,mean(avg_success(find(lickbs_h & success_h),:),1),std(avg_success(find(lickbs_h & success_h),:),[],1)./sqrt(sum(lickbs_h & success_h,2)),'b')
% %     
% %     errorbar(tt,mean(avg_success(find(~lickbs_h & success_h),:),1),std(avg_success(find(~lickbs_h & success_h),:),[],1)./sqrt(sum(~lickbs_h & success_h,2)),'k')
% %     supertitle([date '' mouse ' lick and success: blue (n= ' num2str(sum(lickbs_h & success_h)) ')' ' success but not lick: black (n= ' num2str(sum(~lickbs_h & success_h)) ')']);
% %     saveas(fig, [dest_sub '_lickVSNonLick_Success_avgTCs.fig']);
% %     print([dest_sub '_lickVSNonLick_Success_avgTCs.eps'], '-depsc');
% %     print([dest_sub '_lickVSNonLick_Success_avgTCs.pdf'], '-dpdf');
% %     
% %     fig=figure;
% %     avg_all = [avg_fail avg_success];
% %     ymax = max(max(avg_all,[],2),[],1);
% %     ymin = min(min(avg_all,[],2),[],1);
% %     
% %     if length(notlick_cells) <= 4
% %         ns = 2;
% %     else
% %         ns = ceil(sqrt(length(notlick_cells)));
% %     end
% %     for ic = 1:length(notlick_cells)
% %         n_ic = notlick_cells(ic);
% %         subplot(ns,ceil(length(notlick_cells)/ns),ic)
% %         
% %         hold on;
% %         errorbar(tt,avg_success(n_ic,:), sem_success(n_ic,:),'k')
% %         errorbar(tt,avg_fail(n_ic,:), sem_fail(n_ic,:),'r')
% %         ylim([ymin*1.1 ymax*1.1])
% %         xlim([tt(1) tt(end)])
% %     end
% %     supertitle([date ' ' mouse ' Average success: black (n = ' num2str(size(success_movie,1)) ' trials);' ' fail: red (n= ' num2str(size(fail_movie,1)), ...
% %         ' trials) not lick nCells = ' num2str(length(notlick_cells))])
% %     orient landscape
% %     saveas(fig, [dest_sub '_NotLick_SuccessFail_allTCs.fig']);
% %     print([dest_sub '_NotLick_SuccessFail_allTCs.eps'], '-depsc');
% %     print([dest_sub '_NotLick_SuccessFail_allTCs.pdf'], '-dpdf');
% %     
% %     fig = figure;
% %     hold on
% %     errorbar(tt,mean(avg_fail(find(~lickbs_h & release_h),:),1),std(avg_fail(find(~lickbs_h & release_h),:),[],1)./sqrt(sum(~lickbs_h & release_h,2)),'r')
% %     
% %     errorbar(tt,mean(avg_success(find(~lickbs_h & release_h),:),1),std(avg_success(find(~lickbs_h & release_h),:),[],1)./sqrt(sum(~lickbs_h & release_h,2)),'k')
% %     supertitle([date '' mouse ' NotLick and success: black '  ' NotLick and fail: red (not lick cell n= ' num2str(sum(~lickbs_h & release_h)) ')']);
% %     saveas(fig, [dest_sub '_NotLick_SuccessFail_avgTCs.fig']);
% %     print([dest_sub '_NotLick_SuccessFail_avgTCs.eps'], '-depsc');
% %     print([dest_sub '_NotLick_SuccessFail_avgTCs.pdf'], '-dpdf');
% %     
% %     fig = figure;
% %     success_lick_corr_sig = success_lick_corr(success_lick_corr(:,2) < 0.05, 1);
% %     success_lick_corr_notsig = success_lick_corr(success_lick_corr(:,2) >= 0.05, 1);
% %     fail_lick_corr_sig = fail_lick_corr(fail_lick_corr(:,2) < 0.05, 1);
% %     fail_lick_corr_notsig = fail_lick_corr(fail_lick_corr(:,2) >= 0.05, 1);
% %     hold on
% %     scatter(20*ones(size(success_lick_corr_sig,1),1), success_lick_corr_sig, 20, 'MarkerFaceColor', [0, 0, 0],'MarkerEdgeColor', [0 0 0]);
% %     scatter(20*ones(size(success_lick_corr_notsig,1),1), success_lick_corr_notsig, 20, 'MarkerEdgeColor', [0.5 0.5 0.5]);
% %     scatter(40*ones(size(fail_lick_corr_sig,1),1), fail_lick_corr_sig, 20, 'MarkerFaceColor', [0, 0, 0], 'MarkerEdgeColor', [0 0 0]);
% %     scatter(40*ones(size(fail_lick_corr_notsig,1),1), fail_lick_corr_notsig, 20, 'MarkerEdgeColor', [0.5 0.5 0.5]);
% %     xlim([0 60])
% %     set(gca, 'XTick', [20, 40], 'XTickLabel', {'correct', 'early'});
% %     ylabel('correlation coeffients')
% %     supertitle([date '' mouse ' black: significant correlation found']);
% %     saveas(fig, [dest_sub '_Correlation_lick_dFF_scatter.fig']);
% %     print([dest_sub '_Correlation_lick_dFF_scatter.eps'], '-depsc');
% %     print([dest_sub '_Correlation_lick_dFF_scatter.pdf'], '-dpdf');
% %     
% %     %     fig=figure;
% % %     avg_all = [avg_success];
% % %     ymax = max(max(avg_all,[],2),[],1);
% % %     ymin = min(min(avg_all,[],2),[],1);
% % %     
% % %     for ic = 1:length(success_notlick_cells)
% % %         n_ic = success_notlick_cells(ic);
% % %         n_ic2 = lick_success_cells(ic);
% % %         subplot(2,round(length(success_notlick_cells)/2),ic)
% % %         errorbar(tt,avg_success(n_ic,:), sem_success(n_ic,:),'color', [0.5 0.5 0.5])
% % %         hold on;
% % %         errorbar(tt,avg_success(n_ic2,:), sem_success(n_ic2,:),'k')
% % %         ylim([ymin*1.1 ymax*1.1])
% % %         xlim([tt(1) tt(end)])
% % %     end
% % %     supertitle([date ' ' mouse ' Average lick success: black (n = ' num2str(size(single_lick_movie,1)) ' trials); success not lick resp- gray (n = ' num2str(size(success_movie,...
% % %         1)) ' trials) nCells = ' num2str(length(success_notlick_cells))])
% % %     orient landscape
% % %     saveas(fig, [dest_sub '_SuccessLickcellvsNonlick_allTCs.fig']);
% % %     print([dest_sub '_SuccessLickcellvsNonlick_allTCs.eps'], '-depsc');
% % %     print([dest_sub '_SuccessLickcellvsNonlick_allTCs.pdf'], '-dpdf');
% % end
% % 
% % fig=figure;
% % all_resp = [success_resp_avg fail_resp_avg press_resp_avg];
% % cmax = max(all_resp,[],2);
% % cmin = min(all_resp,[],2);
% % subplot(2,2,1)
% % mask_final_temp = zeros(size(mask_final));
% % mask_final_temp(find(mask_final>0)) = 1;
% % imagesq(reshape(mask_final_temp, [sz(1) sz(2)]));
% % clim([cmin cmax]);
% % set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
% % colorbar
% % 
% % success_mask = mask_final;
% % for ic = 1:nCells
% %     if success_h(ic)
% %         success_mask(find(success_mask==ic))=success_resp_avg(1,ic);
% %     else
% %         success_mask(find(success_mask==ic))=0;
% %     end
% % end
% % success_mask = reshape(success_mask,[sz(1) sz(2)]);
% % subplot(2,2,2)
% % imagesq(success_mask);
% % clim([cmin cmax]);
% % set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
% % title('success')
% % colorbar
% % 
% % fail_mask = mask_final;
% % for ic = 1:nCells
% %     if fail_h(ic)
% %         fail_mask(find(fail_mask==ic))=fail_resp_avg(1,ic);
% %     else
% %         fail_mask(find(fail_mask==ic))=0;
% %     end
% % end
% % fail_mask = reshape(fail_mask,[sz(1) sz(2)]);
% % subplot(2,2,3)
% % imagesq(fail_mask);
% % clim([cmin cmax])
% % set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
% % title('failure')
% % colorbar
% % 
% % press_mask = mask_final;
% % for ic = 1:nCells
% %     if press_h(ic)
% %         press_mask(find(press_mask==ic))=press_resp_avg(1,ic);
% %     else
% %         press_mask(find(press_mask==ic))=0;
% %     end
% % end
% % press_mask = reshape(press_mask,[sz(1) sz(2)]);
% % subplot(2,2,4)
% % imagesq(press_mask);
% % clim([cmin cmax]);
% % set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
% % colorbar
% % title('press')
% % 
% % supertitle([mouse ' ' date ' cell responses']);
% % saveas(fig, [dest_sub '_cell_responses_FOV.fig']);
% % print([dest_sub '_cell_responses_FOV.eps'], '-depsc');
% % print([dest_sub '_cell_responses_FOV.pdf'], '-dpdf');
% % 
% % %plot amplitude response by cell
% % lever = strvcat('press', 'release','success', 'fail');
% % fig=figure;
% % subplot(2,2,1)
% % for ic = 1:nCells
% %     plot(1:4, [press_resp_avg(1,ic) release_resp_avg(1,ic) success_resp_avg(1,ic) fail_resp_avg(1,ic)], '-ok')
% %     hold on
% %     if press_h(1,ic)
% %         plot(1,press_resp_avg(1,ic),'or')
% %         hold on
% %     end
% %     if release_h(1,ic)
% %         plot(2,release_resp_avg(1,ic),'or')
% %         hold on
% %     end
% %     if success_h(1,ic)
% %         plot(3,success_resp_avg(1,ic),'or')
% %         hold on
% %     end
% %     if fail_h(1,ic)
% %         plot(4,fail_resp_avg(1,ic),'or')
% %         hold on
% %     end
% % end
% % set(gca, 'XTick', 1:4, 'XTickLabel', lever);
% % xlim([0.5 4.5])
% % title('Resp amp- red: p<0.05')
% % ylabel('dF/F')
% % 
% % subplot(2,2,2)
% % errorbar(1:4, [mean(press_resp_avg,2) mean(release_resp_avg,2) mean(success_resp_avg,2) mean(fail_resp_avg,2)], [std(press_resp_avg,[],2) std(release_resp_avg,[],2) std(success_resp_avg,[],2) std(fail_resp_avg,[],2)]./sqrt(nCells),'ok')
% % set(gca, 'XTick', 1:4, 'XTickLabel', lever);
% % xlim([0.5 4.5])
% % ylabel('dF/F')
% % title(['Avg resp- all cells: n = ' num2str(nCells)])
% % subplot(2,2,3)
% % errorbar(1:4, [mean(press_resp_avg(1,find(release_h)),2) mean(release_resp_avg(1,find(release_h)),2) mean(success_resp_avg(1,find(release_h)),2) mean(fail_resp_avg(1,find(release_h)),2)], [std(press_resp_avg(1,find(release_h)),[],2) std(release_resp_avg(1,find(release_h)),[],2) std(success_resp_avg(1,find(release_h)),[],2) std(fail_resp_avg(1,find(release_h)),[],2)]./sqrt(sum(release_h,2)),'ok')
% % set(gca, 'XTick', 1:4, 'XTickLabel', lever);
% % xlim([0.5 4.5])
% % ylabel('dF/F')
% % title(['Avg resp- release resp cells: n = ' num2str(sum(release_h,2))])
% % subplot(2,2,4)
% % errorbar(1:4, [mean(press_resp_avg(1,find(press_h)),2) mean(release_resp_avg(1,find(press_h)),2) mean(success_resp_avg(1,find(press_h)),2) mean(fail_resp_avg(1,find(press_h)),2)], [std(press_resp_avg(1,find(press_h)),[],2) std(release_resp_avg(1,find(press_h)),[],2) std(success_resp_avg(1,find(press_h)),[],2) std(fail_resp_avg(1,find(press_h)),[],2)]./sqrt(sum(press_h,2)),'ok')
% % set(gca, 'XTick', 1:4, 'XTickLabel', lever);
% % xlim([0.5 4.5])
% % ylabel('dF/F')
% % title(['Avg resp- press resp cells: n = ' num2str(sum(press_h,2))])
% % supertitle([date ' ' mouse 'avg resp amplitude'])
% % saveas(fig, [dest_sub '_cell_resp_amplitude.fig']);
% % print([dest_sub '_cell_resp_amplitude.eps'], '-depsc');
% % print([dest_sub '_cell_resp_amplitude.pdf'], '-dpdf');
% 
% 
% % success_amp = success_resp - success_base;
% % fail_amp = fail_resp - fail_base;
% % figure;
% % for i = 1:nCells
% %     hold on
% %     scatter(success_amp(:,i), lick_data.lickRateS, 'r');
% %     scatter(fail_amp(:,i), lick_data.lickRateF, 'k');
% % end
% 
