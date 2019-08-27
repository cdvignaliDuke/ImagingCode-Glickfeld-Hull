function [NR_h, NR_p, NR_resp_cells, NR_resp_avg, NR_resp_sem, NR_base, NR_resp] = findRespCell_ITI(NR_movie, pre_cue_frames, ifi)
%uses a ttest to determine significant between a baseline window and a
%response window      NR_movie dim1=trial# dim2=cell# dim3=frame#

%defines baseline and response windows
base_NR_window = pre_cue_frames-floor(750./double(ifi)): pre_cue_frames-round(450./double(ifi));  %lick-500 to lick-200   
resp_NR_window = pre_cue_frames-round(100./double(ifi)) : pre_cue_frames+round(200./double(ifi));  %lick-100 to lick+200

%Averages the values within each window   dim1=trial# dim2=cell#
NR_resp = mean(NR_movie(:,:,resp_NR_window),3);
NR_base = mean(NR_movie(:,:,base_NR_window),3);

NR_resp_avg = mean((NR_resp-NR_base),1);
NR_resp_sem = std((NR_resp-NR_base),[],1)./sqrt(size(NR_resp,1));

%% 2. ttest for significant responses
[NR_h, NR_p] = ttest(NR_resp-NR_base, 0, 'dim', 1, 'tail', 'right');

%% 3. define cells by response to event
NR_resp_cells = find(NR_h);

end