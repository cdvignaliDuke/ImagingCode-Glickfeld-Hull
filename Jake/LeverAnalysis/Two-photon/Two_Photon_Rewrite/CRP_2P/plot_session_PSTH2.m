% plot_session_PSTH2
% script to be run within main2_CRP
% plots PSTHs for each session with sepcific trials removed 
% excluded trials defined in rep_cell_criteria
% plot the session average across neurons and for each neuron

%% load data
%pdf_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\Summary PDF\time courses\';
pdf_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\Summary PDF\PSTH_ex_trials\';
out_dir_TS = ['Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\resp cell troubleshoot\', [session_date, '_', mouse_ID]];
load([dest 'spike_outputs\_evoked_events.mat']);   %loads    3D PSTHs with all cells all trials
load([dest, '\', '_exclude_trials.mat']); %load which trials to exclude 
load([dest, '\', '_pre_cue_ex_cell.mat']);
load([dest '\_cell_categories.mat']);
ifi = 33;
if ~exist(pdf_dir)
    mkdir(pdf_dir);
end
%% assign some variables n shit. 
nCells = length(normalCue);
n = ceil(sqrt(nCells));
if nCells <((n.^2)-n) %function
    n2= n-1;
else
    n2 = n;
end
rew_time = 600;
if session_subset ==1 
    sess_type = 'Day 1 ';
elseif session_subset ==2
    sess_type = 'Post Learning ';
elseif session_subset ==3
    sess_type = 'UR ';
elseif session_subset ==4
    sess_type = '1000ms delay ';
    rew_time = 1100;
end

%% convert to matrices
NR_mat = []; OR_mat = [];  UR_mat = [];
for cell_num = 1:nCells
    NR_mat(cell_num, :,:) = normalCue(cell_num).hist; %dim1=cell dim2=frame dime3=trial#
    if ~isempty(omitCue)
         OR_mat(cell_num, :,:) = omitCue(cell_num).hist; %dim1=cell dim2=frame dime3=trial#
    end
    if ~isempty(unexpCue)
         UR_mat(cell_num, :,:) =unexpCue(cell_num).hist; %dim1=cell dim2=frame dime3=trial#
    end
end
NR_mat = permute(NR_mat, [3,1,2]);
OR_mat = permute(OR_mat, [3,1,2]);
UR_mat = permute(UR_mat, [3,1,2]);

%% exclude specific trials based on criteria from resp_cell_criteria
NR_mean = NaN(nCells, size(NR_mat,3));  %dim1=cell dim2=frame
OR_mean = NaN(nCells, size(NR_mat,3));
UR_mean = NaN(nCells, size(NR_mat,3));
NR_sem = [];
OR_sem = [];
UR_sem = [];
for cell_num = 1:nCells
    temp = squeeze(NR_mat(:,cell_num,:))*30;
    temp(   [NR_ex_trials{cell_num}']  ,:)  = [];
    NR_mean(cell_num,:) = mean(temp);
    NR_sem(cell_num,:) = std(temp)/sqrt(size(temp,1));
    if ~isempty(OR_mat)
        temp = squeeze(OR_mat(:,cell_num,:))*30;
        temp(    [OR_ex_trials{cell_num}']    ,:) = [];
        OR_mean(cell_num,:) = mean(temp);
        OR_sem(cell_num,:) = std(temp)/sqrt(size(temp,1));
    else 
        OR_mean = [];
    end
    if ~isempty(UR_mat)
        temp = squeeze(UR_mat(:,cell_num,:))*30;
        temp( [UR_ex_trials{cell_num}']  ,:) = [];
        UR_mean(cell_num,:) = mean(temp);
        UR_sem(cell_num,:) = std(temp)/sqrt(size(temp,1));
    else 
        UR_mean = [];
    end
end

%%  Plot timecourse for each cell 
x_axis =((-pre_buffer_cue:post_buffer_cue).*double(ifi));
fig=figure('rend', 'painters','pos', [50 100 1000 700]);
avg_all = [NR_mean OR_mean UR_mean];
ymax = max(max(avg_all,[],2),[],1);
ymin = 0;
if ymax <0.1
    ymax = 0.1;
end
%only include resp cells 
resp_cells = find([~pre_cue_ex_cell & allresp_cells_pos]);
for ic = 1:nCells
    if ~ismember(ic, resp_cells);
        continue
    end
    subplot_tight(n,n2,ic);
    if ~isempty(UR_mat)
        errorbar(x_axis, UR_mean(ic,:), UR_sem(ic,:),'g'); hold on;
    end
    if ~isempty(OR_mat)
        errorbar(x_axis, OR_mean(ic,:), OR_sem(ic,:),'r'); hold on;
    end
    errorbar(x_axis, NR_mean(ic,:), NR_sem(ic,:),'k');
    ylim([ymin*1.1 ymax*1.1]);
    xlim([-700 2000]); 
    vline(0,'c'); 
    vline(rew_time, 'b');
end
supertitle([sess_type, mouse_ID, ' PSTH NR=black, OR=red, UR=green. All cells (n=', num2str(nCells), ')']);
orient landscape
print([pdf_dir, sess_type, '_', mouse_ID '_all_PSTHs.pdf'], '-dpdf');

%plot mean TCs for each session
avg_all = [mean(NR_mean,1) mean(OR_mean,1) mean(UR_mean,1)];
ymax = max(max(avg_all,[],2),[],1);
ymin = 0;
if ymax <0.1
    ymax = 0.1;
end
fig2 = figure('rend', 'painters','pos', [50 100 700 700]);
if ~isempty(UR_mean)
    sem_UR_avg = std(UR_mean, [],1)/sqrt(nCells);
    errorbar(x_axis, mean(UR_mean,1), sem_UR_avg,'g'); hold on;
end
if ~isempty(OR_mean)
    sem_OR_avg = std(OR_mean, [],1)/sqrt(nCells);
    errorbar(x_axis, mean(OR_mean,1), sem_OR_avg,'r'); hold on;
end
sem_NR_avg = std(NR_mean, [],1)/sqrt(nCells);
errorbar(x_axis, mean(NR_mean,1), sem_NR_avg,'k'); 
xlim([-700 2000]); 
vline(0, 'c'); 
vline(rew_time, 'b');
ylim([ymin*1.1 ymax*1.1]);
title([sess_type, mouse_ID, ' PSTH NR=black, OR=red, UR=green. All cells (n=', num2str(nCells), ')']);
print([pdf_dir, sess_type, '_', mouse_ID, '_avg_PSTH.pdf'], '-dpdf');


