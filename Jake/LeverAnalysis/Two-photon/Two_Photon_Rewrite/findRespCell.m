function [NR_h, NR_p, NR_resp_cells, NR_resp_avg, NR_resp_sem] = findRespCell(NR_movie, pre_cue_frames, ifi, dest)
%uses a ttest to determine significant between a baseline window and a
%response window

%defines baseline and response windows
base_NR_window = 1:pre_cue_frames-round(200./double(ifi));
resp_NR_window = pre_cue_frames-round(100./double(ifi)):pre_cue_frames+round(2000./double(ifi));

%Averages the values within each window across trials. Then for each cell's
%avg TC it finds the max in the resp window and the min in the baseline
%window. 
[maxx, indp] = max(mean(NR_movie(:,:,resp_NR_window),1),[],3);
[maxx, indc] = min(mean(NR_movie(:,:,base_NR_window),1),[],3);

%define number of cells and allocate memory
nCells = size(NR_movie, 2);
NR_base = zeros(size(NR_movie,1),nCells);
NR_resp = zeros(size(NR_movie,1),nCells);

%main forloop to...
for ic = 1:nCells
    ic
    if indc(ic) - round(50./double(ifi)) <= 0   %indc is the indeces of the min TC values in the baseline window for each cell's avg TC. If the min value is less than 50ms into the TC...
        NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_NR_window(indc(ic)):base_NR_window(indc(ic))+round(100./double(ifi))),3)); %then this baseline value = the avg of the min value and the following 3 frames
    else
     NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_NR_window(indc(ic))-round(50./double(ifi)):base_NR_window(indc(ic))+round(50./double(ifi))),3)); %otherwise the baseline value is determine by a 100ms window surrounding the min val
    end
%     NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_NR_window),3));
    if resp_NR_window(indp(ic))+round(50./double(ifi)) > size(resp_NR_window,2) %if the peak occurs within 50ms of the end of the resp window...
        NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,resp_NR_window(indp(ic))-round(100./double(ifi)):resp_NR_window(indp(ic))),3)); %then the window is taken from the 100ms preceding the peak
    else
        NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,resp_NR_window(indp(ic))-round(50./double(ifi)):resp_NR_window(indp(ic))+round(50./double(ifi))),3)); %otherwise the window is the 100ms surrounding the peak
        
    end
end
NR_resp_avg = mean((NR_resp-NR_base),1);
NR_resp_sem = std((NR_resp-NR_base),[],1)./sqrt(size(NR_resp,1));

save([dest '_cell_resp.mat'], 'NR_base', 'NR_resp');

%% 2. ttest for significant responses

[NR_h, NR_p] = ttest(NR_base, NR_resp, 'dim', 1, 'tail', 'both');

%% 3. define cells by response to event

NR_resp_cells = find(NR_h);
end