%Plot TCs and convert to PDFs
%Store in PDF sumamry folder. 
%Script to be run with main_CRP2

%% load data
%pdf_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\Summary PDF\time courses\';
pdf_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\Summary PDF\time courses_ex_trials\';
%out_dir_TS = ['Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\resp cell troubleshoot\', [session_date, '_', mouse_ID]];
%load([dest '_cell_TCs.mat']);  %2D mats of TCs averaged across trials
load([dest '_cue_movies.mat']);   %loads NR_movie etc   3D TCs with all cells all trials
load([dest, '\', '_exclude_trials.mat']); %load which trials to exclude 
load([dest '\_cell_categories.mat']);
load([dest, '\', '_pre_cue_ex_cell.mat']); %load which cells to exclude
if ~exist('ex_cell_all')
    ex_cell_all = [];
end
temp_cell = [sum(pre_cue_ex_cell); length(pre_cue_ex_cell)];
if ~exist(pdf_dir)
    mkdir(pdf_dir);
end

%% assign some variables n shit. 
nCells = size(NR_movie,2);
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

%% exclude specific trials based on criteria from resp_cell_criteria
NR_mean = NaN(nCells, size(NR_movie,3));  %dim1=cell dim2=frame
OR_mean = NaN(nCells, size(NR_movie,3));
UR_mean = NaN(nCells, size(NR_movie,3));
NR_sem = [];
OR_sem = [];
UR_sem = [];
for cell_num = 1:nCells
    temp = squeeze(NR_movie(:,cell_num,:));
    temp(    unique(sort([NR_ex_trials{cell_num}', NR_ex_trials{cell_num}']))   ,:)  = [];
    NR_mean(cell_num,:) = mean(temp);
    NR_sem(cell_num,:) = std(temp)/sqrt(size(temp,1));
    if ~isempty(OR_movie)
        temp = squeeze(OR_movie(:,cell_num,:));
        temp(    unique(sort([OR_ex_trials{cell_num}', OR_ex_trials{cell_num}']))    ,:) = [];
        OR_mean(cell_num,:) = mean(temp);
        OR_sem(cell_num,:) = std(temp)/sqrt(size(temp,1));
    else 
        OR_mean = [];
    end
    if ~isempty(UR_movie)
        temp = squeeze(UR_movie(:,cell_num,:));
        temp( unique(sort([UR_ex_trials{cell_num}', UR_ex_trials{cell_num}']))  ,:) = [];
        UR_mean(cell_num,:) = mean(temp);
        UR_sem(cell_num,:) = std(temp)/sqrt(size(temp,1));
    else 
        UR_mean = [];
    end
end

%%  Plot timecourse for each cell 
x_axis =((-pre_cue_frames:post_cue_frames).*double(ifi));
fig=figure('rend', 'painters','pos', [50 100 1000 700]);
avg_all = [NR_mean OR_mean UR_mean];
ymax = max(max(avg_all,[],2),[],1);
ymin = min(min(avg_all,[],2),[],1);
if ymax <0.1
    ymax = 0.1;
end
if ymin >-0.1
    ymin = -0.1;
end
for ic = 1:nCells
    subplot_tight(n,n2,ic);
    if ismember(ic, find(pre_cue_ex_cell))
        continue
    end
    if ~ismember(ic, find(rew_cells_pos)) %NR_Rew_resp_cells_pos
        continue
    end
    if ~isempty(UR_movie)
        errorbar(x_axis, UR_mean(ic,:), UR_sem(ic,:),'g'); hold on;
    end
    if ~isempty(OR_movie)
        errorbar(x_axis, OR_mean(ic,:), OR_sem(ic,:),'r'); hold on;
    end
    errorbar(x_axis, NR_mean(ic,:), NR_sem(ic,:),'k');
    ylim([ymin*1.1 ymax*1.1]);
    xlim([-700 2000]); 
    vline(0,'c'); 
    vline(rew_time, 'b');
end
supertitle([sess_type, mouse_ID, ' timecourse NR=black, OR=red, UR=green. All cells (n=', num2str(nCells), ') (non NR rew +resp)']);
orient landscape
print([pdf_dir, sess_type, '_', mouse_ID '_all_TCs_non NR rew +resp.pdf'], '-dpdf');

%plot mean TCs for each session
avg_all = [mean(NR_mean,1) mean(NR_mean,1) mean(NR_mean,1)];
ymax = max(max(avg_all,[],2),[],1);
ymin = min(min(avg_all,[],2),[],1);
if ymax <0.1
    ymax = 0.1;
end
if ymin >-0.1
    ymin = -0.1;
end 

fig2 = figure('rend', 'painters','pos', [50 100 700 700]);
for plot_num = 1:2
    if plot_num == 2
        cell_use_logical = (~pre_cue_ex_cell) & rew_cells_pos; %
        subplot(1,2,2);
    else
        cell_use_logical = ones(1,nCells);
        subplot(1,2,1);
    end
    if ~isempty(UR_mean)
        UR_mean = UR_mean(logical([cell_use_logical]),:);  
        sem_UR_avg = std(UR_mean, [],1)/sqrt(nCells);
        errorbar(x_axis, mean(UR_mean,1), sem_UR_avg,'g'); hold on;
    end
    if ~isempty(OR_mean)
        OR_mean = OR_mean(logical([cell_use_logical]),:);
        sem_OR_avg = std(OR_mean, [],1)/sqrt(nCells);
        errorbar(x_axis, mean(OR_mean,1), sem_OR_avg,'r'); hold on;
    end
    NR_mean = NR_mean(logical([cell_use_logical]),:);
    sem_NR_avg = std(NR_mean, [],1)/sqrt(nCells);
    errorbar(x_axis, mean(NR_mean,1), sem_NR_avg,'k');
    xlim([-700 2000]);
    vline(0,'c');
    vline(rew_time, 'b');
    ylim([ymin*1.1 ymax*1.1]);
    if plot_num ==1 
        title(['All cells: (n=', num2str(sum(cell_use_logical)), ')']);
    else
        title(['Resp Cells (n=', num2str(sum(cell_use_logical)), ')']);
    end
end
suptitle([sess_type,' ', mouse_ID, ' timecourse NR=black, OR=red, UR=green.  (non NR rew +resp)']);
print([pdf_dir, sess_type, '_', mouse_ID, '_avg_TCs_non NR rew +resp.pdf'], '-dpdf');



