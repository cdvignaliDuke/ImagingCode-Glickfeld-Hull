function [NR_h, NR_p, NR_resp_cells, NR_resp_avg, NR_resp_sem, NR_base, NR_resp] = findRespCell(NR_movie, event_frame, ifi, base_buffer, resp_buffer, effect_sign)
%function for identifying significantly responsive neurons using the df/f traces from dendritic masks
%NR_movie   dim1=trials  dim2=cells  dim3=frames

if isempty(NR_movie)
    NR_h =[]; NR_p =[]; NR_resp_cells =[]; NR_resp_avg =[]; NR_resp_sem =[]; NR_base =[]; NR_resp =[];
    return
end

%define baseline and response windows
base_NR_window = event_frame-round(base_buffer./double(ifi)):event_frame -1;
resp_NR_window = event_frame:event_frame+round(resp_buffer./double(ifi));

%find the min and max value indeces in the base and resp windows respectively
[~, indp] = max(nanmean(NR_movie(:,:,resp_NR_window),1),[],3);
[~, indp_neg] = min(nanmean(NR_movie(:,:,resp_NR_window),1),[],3);
    
%replace with code that takes the mean of whole baseline window
[~, indc] = min(nanmean(NR_movie(:,:,base_NR_window),1),[],3);

%allocate memory and define variables
nCells = size(NR_movie, 2);
NR_base = zeros(size(NR_movie,1),nCells);
NR_resp = zeros(size(NR_movie,1),nCells);
NR_resp_neg = zeros(size(NR_movie,1),nCells);

%for each cell...
for ic = 1:nCells
    %find mean baseline response
    NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,[base_NR_window]),3));  

    %define response window and store value  POSITIVE
    if resp_NR_window(indp(ic))+round(50./double(ifi)) > size(resp_NR_window,2) %make sure there is enough room in the rest of the TC for the window. 
        NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,resp_NR_window(indp(ic))-round(100./double(ifi)):resp_NR_window(indp(ic))),3));
    else
        NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,resp_NR_window(indp(ic))-round(50./double(ifi)):resp_NR_window(indp(ic))+round(50./double(ifi))),3));
    end
    
    %define response window and store value  NEGATIVE
    if resp_NR_window(indp_neg(ic))+round(50./double(ifi)) > size(resp_NR_window,2) %make sure there is enough room in the rest of the TC for the window. 
        NR_resp_neg(:,ic) = squeeze(mean(NR_movie(:,ic,resp_NR_window(indp_neg(ic))-round(100./double(ifi)):resp_NR_window(indp_neg(ic))),3));
    else
        NR_resp_neg(:,ic) = squeeze(mean(NR_movie(:,ic,resp_NR_window(indp_neg(ic))-round(50./double(ifi)):resp_NR_window(indp_neg(ic))+round(50./double(ifi))),3));
    end
end

%================================================
NR_resp_avg = []; %nanmean((NR_resp-NR_base),1);
NR_resp_sem = [];  %nanstd((NR_resp-NR_base),[],1)./sqrt(size(NR_resp,1));

%% 2. ttest for significant responses
if isempty(find(~isnan(NR_base)))
    NR_h = zeros(1,size(NR_base,2));
    NR_p = zeros(1,size(NR_base,2));
else
    if strcmp(effect_sign, 'pos') ==1
        [NR_h, NR_p] = ttest(NR_resp, NR_base, 'dim', 1, 'tail', 'right', 'alpha', 0.05);
    elseif strcmp(effect_sign, 'neg') ==1
        [NR_h, NR_p] = ttest(NR_resp_neg, NR_base, 'dim', 1, 'tail', 'left', 'alpha', 0.05);
    end
end
%% 3. define cells by response to event
NR_h = logical(NR_h);
NR_resp_cells = find(NR_h);

end


