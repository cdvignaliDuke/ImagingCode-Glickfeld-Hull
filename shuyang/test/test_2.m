% use other criteria to sort out bad cells:
%e.g. session:191115_img1041
peak_std = zeros(1,size(spk_peak,2));
peak_mean = zeros(1,size(spk_peak,2));
for c = 1:size(spk_peak,2)                  % for each cell
    peak_std(c) = std(spk_peak{c});
    peak_mean(c) = mean(spk_peak{c});
end
peak_variation = peak_std./peak_mean; %if peak variation is too big, this cell might not be good.
figure; plot(peak_variation); 
title('peak variation');
peak_vartoomuch = find(peak_variation>0.64);
% cannot totally through out all of the cells in here, need to go back and
% look at each cell that has a high variation and decide.
% but this at least gives me a good pool to look at
% cell # 18,9,7,4,1 can all be sorted out using this way
% ----------------------------------------------------------------------------------------------
% another feature is that the deconvolved signal doesn't go back to
% baseline for a long time: cell#88,53,91,and 101
nobase_maxperiod = zeros(1,size(spk_peak,2));
for c = 1:size(spk_peak,2)
    nobase_maxperiod(c) = max(diff(find(kernel(:,c)<100))); %kernel values smaller than 100 is considered as basline. diff(find) gives you how many continous frames the denoised signal doesn't return to baseline
end
figure; plot(nobase_maxperiod);
title('not return to baseline period');
%if the maximum not return to baseline period is bigger than n # of frames,
%through out this cell
nobase_thres = 1000;
cells = 1:size(spk_peak,2);
nobase_cell = find(nobase_maxperiod > nobase_thres);

%----------------------------------------------------------------------------------------------------
%a third feature is that the cells are firing a lot in a short period of
%time
%calculate the firing rate of a period of time and if the FR is high for a
%long time then through the cell out first calculate # of spikes per second
maxfire = zeros(1,size(spk_peak,2));
binwidth = 900; % number of frames
bins = floor(size(spk_logic,1)/binwidth);
for c = 1:size(spk_peak,2)
    nfire = []; t = 1;
    while t < size(spk_logic,1)-binwidth
        nfire = [nfire sum(spk_logic(t:t+binwidth-1,c)==1)];
        t = t+1;
    end
    maxfire(c) = max(nfire);
end
figure; plot(maxfire); title('firing too much');
figure; hist(maxfire); 
FR_toohigh = find(maxfire>65);

badcells2 = union(peak_vartoomuch,nobase_cell);
badcells2 = union(badcells2,FR_toohigh);










