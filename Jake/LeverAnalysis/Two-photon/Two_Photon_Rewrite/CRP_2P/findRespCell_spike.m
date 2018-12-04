function [NR_h, NR_p, NR_resp_cells, NR_resp_avg, NR_resp_sem, NR_base, NR_resp] = findRespCell_spike(NR_movie, base_NR_window, resp_NR_window, effect_sign)
%function for identifying significantly responsive neurons using the Ca event PSTHs from dendritic masks
%NR_movie  is a structure with fields for multiple variables for each cell including hist (histograms)

%output empty matrices if this session does not have any trials of this trial type
if isempty(NR_movie) 
    NR_h =[]; NR_p =[]; NR_resp_cells =[]; NR_resp_avg =[]; NR_resp_sem =[]; NR_base =[]; NR_resp =[];
    return
end

%extract histogram data into a matrix
NR_movie_hist = cat(3, NR_movie.hist); %dim1=frame#  dim2=trial#  dim3=cell#

%allocate memory and define variables
nCells = size(NR_movie_hist, 3);
NR_base = zeros(size(NR_movie_hist,2),nCells);
NR_resp = zeros(size(NR_movie_hist,2),nCells);
NR_resp_neg = zeros(size(NR_movie_hist,2),nCells);

%for each cell...
for ic = 1:nCells
    %find mean baseline 
    NR_base(:,ic) = squeeze(nanmean(NR_movie_hist([base_NR_window],:,ic),1));  
    %find mean response  
    NR_resp(:,ic) = squeeze(nanmean(NR_movie_hist([resp_NR_window],:,ic),1));  
end

%find the mean resp size across trials for each cell
NR_resp_avg = nanmean((NR_resp-NR_base),1);
NR_resp_sem = nanstd((NR_resp-NR_base),[],1)./sqrt(size(NR_resp,1));

%% 2. ttest for significant responses
if isempty(find(~isnan(NR_base)))
    NR_h = zeros(1,size(NR_base,2));
    NR_p = zeros(1,size(NR_base,2));
else
    if strcmp(effect_sign, 'pos') ==1
        [NR_h, NR_p] = ttest(NR_resp, NR_base, 'dim', 1, 'tail', 'right', 'alpha', 0.05);
    elseif strcmp(effect_sign, 'neg') ==1
        [NR_h, NR_p] = ttest(NR_resp, NR_base, 'dim', 1, 'tail', 'left', 'alpha', 0.05);
    end
end

%% 3. define cells by response to event
NR_h = logical(NR_h);
NR_resp_cells = find(NR_h);

end


