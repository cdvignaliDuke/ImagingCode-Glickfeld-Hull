function [movie_mean, movie_sem, movie_mean_all, movie_sem_all] = get_movie_mean_sem(R_movie, trial_cond, lick_cond);
%function for finding the mean and sem first across trials then across
%trials and cells for the df/f TCs extracted from the dendritic masks
% IF condition for selecting the first half of trials only
% IF condition to select all trials vs no lick trials

%average and sem across Trials
if strcmp(trial_cond, 'all') | isempty(trial_cond)
    movie_mean = squeeze(nanmean(R_movie,1));
elseif strcmp(trial_cond, 'half')
    nTrial = size(R_movie,1);
    movie_mean = squeeze(nanmean(R_movie(1:round(nTrial/2),:,:),1));
end
if strcmp(lick_cond, 'all') | isempty(lick_cond)
    movie_sem = squeeze(std(R_movie,1)./sqrt(size(R_movie,1)));
elseif strcmp(lick_cond, 'nolick')
    movie_mean_temp = R_movie(~isnan(R_movie(:,1,1)),:,:);
    movie_sem = squeeze(nanstd(R_movie,1)./sqrt(size(movie_mean_temp,1)));
end

%average and sem across cells and trials
movie_mean_all = mean(movie_mean,1);
movie_sem_all = std(movie_mean,1)./sqrt(size(movie_mean,1));


end



