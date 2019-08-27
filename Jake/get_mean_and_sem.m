function [data_mean, data_sem] = get_mean_and_sem(input_data);

if isempty(find(isnan(input_data)))
    data_mean = mean(input_data);
    if size(input_data,1) ~= 1 & size(input_data,2)~= 1
        data_sem = std(input_data)/sqrt(size(input_data,1));
    else
        data_sem = std(input_data)/sqrt(length(input_data));
    end
else
    data_mean = nanmean(input_data);
    if size(input_data,1) ~= 1 & size(input_data,2)~= 1
        data_sem = std(input_data([isfinite(input_data)])) / sqrt(sum(~isnan(input_data),[],1));
    else
        data_sem = std(input_data([isfinite(input_data)])) / sqrt(sum(~isnan(input_data)));
    end
end

   
return
