%TC_diag_criteria
%diagnose trial TCs. Look at single neuron single trial responses. Use these plots to generate 
% and test criteria for selecting which neurons to include/exclude.
% to be run within main_CRP2
%NR_movie dim1=trials  dim2=cell# dim3=frame#

%load data
load([dest '_cue_movies.mat']);
x_axis = [-pre_cue_frames:post_cue_frames]*double(ifi);

%have user select a specific neuron. 
%neuron_id = input('Enter the neuron # to plot: ');
for neuron_id = 1:15
NR_cell_mean = squeeze(NR_movie(:,neuron_id,:));
this_cell_trials = NR_cell_mean;
OR_cell_mean = [];
if ~isempty(OR_movie)
    OR_cell_mean = squeeze(OR_movie(:,neuron_id,:));
    this_cell_trials = [this_cell_trials; OR_cell_mean];
end
UR_cell_mean= [];
if~isempty(UR_movie)
    UR_cell_mean = squeeze(UR_movie(:,neuron_id,:));
    this_cell_trials = [this_cell_trials; UR_cell_mean];
end

%Plot the mean and 3 std of all trials for that neuron
figure; 
subplot(2,1,1);
shadedErrorBar(x_axis, mean(this_cell_trials), std(this_cell_trials,[],1)*4);
hold on;
for trial_num=1:size(NR_cell_mean,1)
    if trial_num<size(NR_cell_mean,1)/2
        plot(x_axis, NR_cell_mean(trial_num,:),'k');
    else
        plot(x_axis, NR_cell_mean(trial_num,:),'Color', [0.1, 0.1, 0.1]);
    end
end
ylim([-0.5 1]);

subplot(2,1,2);
shadedErrorBar(x_axis, mean(this_cell_trials), std(this_cell_trials,[],1)*4);
hold on;
if ~isempty(OR_cell_mean)
    for trial_num=1:size(OR_cell_mean,1)
        if trial_num<size(OR_cell_mean,1)/2
        plot(x_axis, OR_cell_mean(trial_num,:),'r');
        else
            plot(x_axis, OR_cell_mean(trial_num,:), 'Color', [1,0.1,0.6]);
        end
    end
end
ylim([-0.5 1])
title([mouse_num, '-', session_date, '-neuron ', neuron_id]);

% if ~isempty(UR_cell_mean)
%     for trial_num=1:size(UR_cell_mean,1)
%         plot(x_axis, UR_cell_mean(trial_num,:),'k');
%     end
% end
end
%include vlines for cue/rew

%then plot all single trials on top of that. Color coded for NR/OR/UR






