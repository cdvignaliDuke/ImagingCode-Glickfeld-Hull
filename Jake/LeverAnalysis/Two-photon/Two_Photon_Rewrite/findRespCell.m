function [NR_h, NR_p, NR_resp_cells, NR_resp_avg, NR_resp_sem] = findRespCell(NR_movie, pre_cue_frames, ifi, dest)
base_NR_window = 1:pre_cue_frames-round(200./double(ifi));
resp_NR_window = pre_cue_frames-round(100./double(ifi)):pre_cue_frames+round(2000./double(ifi));
[maxx, indp] = max(mean(NR_movie(:,:,resp_NR_window),1),[],3);
[maxx, indc] = min(mean(NR_movie(:,:,base_NR_window),1),[],3);
nCells = size(NR_movie, 2);
NR_base = zeros(size(NR_movie,1),nCells);
NR_resp = zeros(size(NR_movie,1),nCells);
for ic = 1:nCells
    ic
    if indc(ic) - round(50./double(ifi)) <= 0
        NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_NR_window(indc(ic)):base_NR_window(indc(ic))+round(100./double(ifi))),3));
    else
     NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_NR_window(indc(ic))-round(50./double(ifi)):base_NR_window(indc(ic))+round(50./double(ifi))),3));
    end
%     NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_NR_window),3));
    if resp_NR_window(indp(ic))+round(50./double(ifi)) > size(resp_NR_window,2)
        NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,resp_NR_window(indp(ic))-round(100./double(ifi)):resp_NR_window(indp(ic))),3));
    else
        NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,resp_NR_window(indp(ic))-round(50./double(ifi)):resp_NR_window(indp(ic))+round(50./double(ifi))),3));
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