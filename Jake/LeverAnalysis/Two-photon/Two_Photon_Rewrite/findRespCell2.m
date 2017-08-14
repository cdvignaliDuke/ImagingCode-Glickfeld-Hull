function [NR_h, NR_p, NR_resp_cells, NR_resp_avg, NR_resp_sem, NR_base, NR_resp] = findRespCell2(NR_movie, pre_cue_frames, ifi)
%uses a ttest to determine significant between a baseline window and a
%response window

%defines baseline and response windows
base_NR_window = pre_cue_frames-round(500./double(ifi)): pre_cue_frames-round(200./double(ifi));
resp_NR_window = pre_cue_frames-round(100/double(ifi)) : pre_cue_frames+round(100./double(ifi));

%Averages the values within each window across trials. Then compare values
%between a resp window and a baseline window
NR_resp = mean(mean(NR_movie(:,:,resp_NR_window),1),3);
NR_base = mean(mean(NR_movie(:,:,base_NR_window),1),3);

NR_resp_avg = mean((NR_resp-NR_base),1);
NR_resp_sem = std((NR_resp-NR_base),[],1)./sqrt(size(NR_resp,1));


%% 2. ttest for significant responses
[NR_h, NR_p] = ttest(NR_base, NR_resp, 'dim', 1, 'tail', 'both');

%% 3. define cells by response to event
NR_resp_cells = find(NR_h);

end