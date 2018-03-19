clear
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
%days = {'150518_img24', '150519_img24', '150518_img25', '150517_img25', '150716_img27', '150718_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41'}; %'150718_img27', '150719_img27',
%days = {'150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41'}; %'150718_img27', '150719_img27',
%days = {'151009_img30', '151011_img30', '151212_img32', '160129_img36', '160314_img38', '160315_img38', '160606_img46'}; %only the days which meet bx criteria
%days = {'150716_img28', '150717_img28', '151021_img29', '151022_img29', '151211_img32', '160129_img35', '160131_img35', '160131_img36', '160319_img41', '160320_img41'}; %Only the days which do NOT meet bx criteria
WF_plotting_lists_of_days;
%days = {'150716_img28'};
reactTimesMs_all = [];
reqHoldTimeMs_all = [];
across_animals_RTs_binned=[];
bin_duration =250;
%figure;
hold on; subplot(2,1,2);
for kk=1:length(days)
    day_name  =  days{kk};
    bfile = dir([BEHAVE_DIR 'data-*i9' days{kk}(end-1:end) '-' days{kk}(1:6) '*' ]);
    behave_dest = [BEHAVE_DIR bfile.name];
    b_data = load(behave_dest);
    %if b_data.input.randReqHoldMaxMs > 1501;
%     if b_data.input.randReqHoldMaxMs < 4500;
%         continue
%     end
%     if b_data.input.fixedReqHoldTimeMs < 499 %subset of sessions with a 400ms fixed req hold time is throwing off the fit
%         continue
%     end
    reactTimesMs = cell2mat(b_data.input.reactTimesMs);
    reqHoldTimeMs = round(cell2mat(b_data.input.tTotalReqHoldTimeMs));
    %assert(isempty(find(reqHoldTimeMs<499)))
    %only look at reaction times for correct trials
    corr_inx = find(reactTimesMs>150 & reactTimesMs<999);
    reactTimesMs = reactTimesMs(corr_inx);
    reqHoldTimeMs = reqHoldTimeMs(corr_inx);
    %scatter(reqHoldTimeMs, reactTimesMs);
    reactTimesMs_all = [reactTimesMs_all, reactTimesMs];
    reqHoldTimeMs_all = [reqHoldTimeMs_all, reqHoldTimeMs];
    model_1 = fitlm(reqHoldTimeMs,double(reactTimesMs));
    plot(model_1); hold on;
    
    %for each animal, bin the reqHoldDur into 250ms bins then take the mean RT for each bin for that animal.
    curr_edge_min = 0;
    curr_edge_max = bin_duration;
    num_bins = 5500/bin_duration;
    for this_bin = 1:num_bins
        holds_this_bin = find(reqHoldTimeMs>curr_edge_min & reqHoldTimeMs<=curr_edge_max);
        if isempty(holds_this_bin)
            RT_binned(this_bin) = NaN;
        else
            RTs_this_bin = reactTimesMs([holds_this_bin]);
            RT_binned(this_bin) = mean(RTs_this_bin);
        end
        curr_edge_min = curr_edge_max;
        curr_edge_max = curr_edge_max + bin_duration;
    end
    across_animals_RTs_binned(:,end+1) = RT_binned';
end
refline(0,0);
title('reaction times vs total required hold times: long duration trials');
ylabel('Reaction time relative to the visual cue');
xlabel('total required hold time');
model_1 = fitlm(reqHoldTimeMs_all,double(reactTimesMs_all));
figure; 
plot(model_1);


%plot RT binned according to reqHold
across_animals_RTs_binned(2,:) = NaN;
figure; 
x_axis = [1:22]*bin_duration;
%x_axis = [0:bin_duration:5250];
errorbar(x_axis, nanmean(across_animals_RTs_binned,2), nanstd(across_animals_RTs_binned, [], 2)./sqrt(  sum(~isnan(across_animals_RTs_binned),2)   ),'-ok');
ylim([0 450]); xlim([450 5100]);
ylabel('reaction time (ms)');
xlabel('required hold duration (250ms bins)');
title(['Required hold duration vs reaction time for correct trials: n=', num2str(size(across_animals_RTs_binned,2))]);



