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
[~, indc] = min(nanmean(NR_movie(:,:,base_NR_window),1),[],3);

%allocate memory and define variables
nCells = size(NR_movie, 2);
NR_base = zeros(size(NR_movie,1),nCells);
NR_resp = zeros(size(NR_movie,1),nCells);

%for each cell...
for ic = 1:nCells
    %define the baseline window and store value
    if indc(ic) - round(50./double(ifi)) <= 0  %if the min value is 50ms or less from the start of the TC..
        NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_NR_window(indc(ic)):base_NR_window(indc(ic))+round(100./double(ifi))),3)); %use the mean of 100ms after min frame as baseline 
    else
        NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_NR_window(indc(ic))-round(50./double(ifi)):base_NR_window(indc(ic))+round(50./double(ifi))),3)); %else use a window -50:50ms around min frame
    end
    %define response window and store value
    if resp_NR_window(indp(ic))+round(50./double(ifi)) > size(resp_NR_window,2) %make sure there is enough room in the rest of the TC for the window. 
        NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,resp_NR_window(indp(ic))-round(100./double(ifi)):resp_NR_window(indp(ic))),3));
    else
        NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,resp_NR_window(indp(ic))-round(50./double(ifi)):resp_NR_window(indp(ic))+round(50./double(ifi))),3));
    end
end
NR_resp_avg = nanmean((NR_resp-NR_base),1);
NR_resp_sem = nanstd((NR_resp-NR_base),[],1)./sqrt(size(NR_resp,1));

%% 2. ttest for significant responses
if strcmp(effect_sign, 'pos') ==1
    [NR_h, NR_p] = ttest(NR_resp, NR_base, 'dim', 1, 'tail', 'right');
elseif strcmp(effect_sign, 'neg') ==1
    [NR_h, NR_p] = ttest(NR_resp, NR_base, 'dim', 1, 'tail', 'left');
elseif strcmp(effect_sign, 'both') ==1 | isempty(effect_sign)
    [NR_h, NR_p] = ttest(NR_resp, NR_base, 'dim', 1, 'tail', 'both');
end
    
%% 3. define cells by response to event
NR_h = logical(NR_h);
NR_resp_cells = find(NR_h);

end


