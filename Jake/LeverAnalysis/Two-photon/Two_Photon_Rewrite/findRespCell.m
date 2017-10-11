function [NR_h, NR_p, NR_resp_cells, NR_resp_avg, NR_resp_sem, NR_base, NR_resp, NR_pos_h, NR_neg_h] = findRespCell(NR_movie, event_frame, ifi, base_buffer, resp_buffer)
%NR_movie dim1=trial#  dim2=cell#  dim3=frame#
%uses inputs to determine analysis windows for base and response values
base_NR_window = event_frame-round(base_buffer./double(ifi)):event_frame -1;
resp_NR_window = event_frame:event_frame+round(resp_buffer./double(ifi));
%avgs together all the trials then finds the index of the max/min resp for each cell
[~, indp] = max(nanmean(NR_movie(:,:,resp_NR_window),1),[],3);
[~, indc] = min(nanmean(NR_movie(:,:,base_NR_window),1),[],3);
nCells = size(NR_movie, 2);
NR_base = zeros(size(NR_movie,1),nCells);
NR_resp = zeros(size(NR_movie,1),nCells);

for ic = 1:nCells
    if indc(ic) - round(50./double(ifi)) <= 0 %if the min resp for this cell happens within 50ms of the start of the TC...
        NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_NR_window(indc(ic)):base_NR_window(indc(ic))+round(100./double(ifi))),3)); %then use the 100ms following the index value
    else %else use the 50ms before and after the indexed value
     NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_NR_window(indc(ic))-round(50./double(ifi)):base_NR_window(indc(ic))+round(50./double(ifi))),3));
    end %     NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_NR_window),3));
    
    if resp_NR_window(indp(ic))+round(50./double(ifi)) > size(resp_NR_window,2)
        NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,resp_NR_window(indp(ic))-round(100./double(ifi)):resp_NR_window(indp(ic))),3));
    else
        NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,resp_NR_window(indp(ic))-round(50./double(ifi)):resp_NR_window(indp(ic))+round(50./double(ifi))),3));
    end
end
NR_resp_avg = nanmean((NR_resp-NR_base),1);
NR_resp_sem = nanstd((NR_resp-NR_base),[],1)./sqrt(size(NR_resp,1));

% save([dest '_cell_resp.mat'], 'NR_base', 'NR_resp');

%% 2. ttest for significant responses

[NR_h, NR_p] = ttest(NR_base, NR_resp, 'dim', 1, 'tail', 'both');

[NR_pos_h, NR_pos_p] = ttest(NR_resp - NR_base, 0, 'dim', 1, 'tail', 'right');

[NR_neg_h, NR_neg_p] = ttest(NR_resp - NR_base, 0, 'dim', 1, 'tail', 'left');

if ~isempty(NR_movie)
    NR_h = logical(NR_h);
    NR_pos_h = logical(NR_pos_h);
    NR_neg_h = logical(NR_neg_h);
end
%% 3. define cells by response to event

NR_resp_cells = find(NR_h);
end