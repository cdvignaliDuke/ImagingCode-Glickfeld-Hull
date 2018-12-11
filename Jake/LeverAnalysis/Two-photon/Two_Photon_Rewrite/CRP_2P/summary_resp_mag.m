% Scropt for determining the magnitude of cell responses
clear
file_info_CRP_all
data_base = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\';
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\';

%% Day 1
D1_cell_cats=[];  D1_pvals=[];  D1_cell_resp=[]; D1_ex_cell=[];
for sub = 1:size(days_1,2)
    %collect session/mouse information  set pathnames
    if isempty(days_1{sub})
        continue
    end
    [dest_sub, ~, ~, ~, session, session_date, mouse_num, mouse_ID, rID] = get_sess_info(days_1{sub}, runID, crp_dir, []);
    
    %load data into cell
    D1_cell_cats{sub} = load([dest_sub '_cell_categories.mat']);
    D1_cell_resp{sub} = load([dest_sub '_cell_resp.mat']);
    D1_ex_cell{sub} = load([dest_sub, '_pre_cue_ex_cell.mat']);
end

%% Post learning
Dpost_cell_cats = [];  Dpost_pvals = [];  Dpost_cell_resp = []; Dpost_ex_cell=[];
for sub = 1:size(days_post,2)
    %collect session/mouse information  set pathnames
    [dest_sub, ~, ~, ~, session, session_date, mouse_num, mouse_ID, rID] = get_sess_info(days_post{sub}, runID, crp_dir, []);

    %load data into cell
    Dpost_cell_cats{sub} = load([dest_sub '_cell_categories.mat']);
    Dpost_cell_resp{sub} = load([dest_sub '_cell_resp.mat']);
    Dpost_ex_cell{sub} = load([dest_sub, '_pre_cue_ex_cell.mat']);
end

%% Unexpected Reward
UR_cell_cats = [];  UR_pvals = [];  UR_cell_resp = []; UR_ex_cell=[];
for sub = 1:size(days_UR,2)
    %collect session/mouse information  set pathnames
    [dest_sub, ~, ~, ~, session, session_date, mouse_num, mouse_ID, rID] = get_sess_info(days_UR{sub}, runID, crp_dir, []);

    %load data into cell
    UR_cell_cats{sub} = load([dest_sub '_cell_categories.mat']);
    UR_cell_resp{sub} = load([dest_sub '_cell_resp.mat']);
    UR_ex_cell{sub} = load([dest_sub, '_pre_cue_ex_cell.mat']);
end

%% 1000ms delay
D1000_cell_cats=[];  D1000_pvals=[];  D1000_cell_resp=[]; D1000_ex_cell=[];
for sub = 1:size(days_1000,2)
    %collect session/mouse information  set pathnames
    if isempty(days_1000{sub})
        continue
    end
    [dest_sub, ~, ~, ~, session, session_date, mouse_num, mouse_ID, rID] = get_sess_info(days_1000{sub}, runID, crp_dir, []);
    
    %load data into cell
    D1000_cell_cats{sub} = load([dest_sub '_cell_categories.mat']);
    D1000_cell_resp{sub} = load([dest_sub '_cell_resp.mat']);
    D1000_ex_cell{sub} = load([dest_sub, '_pre_cue_ex_cell.mat']);
end

%% Aggregate response magnitude data
D1_cue_resp_mag = [];
D1_rew_resp_mag = [];
for sub = 1:length(D1_cell_resp)
    if isempty(days_1{sub})
        continue
    end
    use_cells = [~D1_ex_cell{sub}.pre_cue_ex_cell & D1_cell_cats{sub}.allresp_cells];
    %use_cells = ~D1_ex_cell{sub}.pre_cue_ex_cell;
    D1_cue_resp_mag = [D1_cue_resp_mag, nanmean(D1_cell_resp{sub}.NR_Cue_resp_pos(:,[use_cells])) - nanmean(D1_cell_resp{sub}.NR_Cue_base_pos(:,[use_cells]))];
    D1_rew_resp_mag = [D1_rew_resp_mag, nanmean(D1_cell_resp{sub}.NR_Rew_resp_pos(:,[use_cells])) - nanmean(D1_cell_resp{sub}.NR_Rew_base_pos(:,[use_cells]))];
end
PL_cue_resp_mag = [];
PL_rew_resp_mag = [];
for sub = 1:length(Dpost_cell_resp)
    use_cells = [~Dpost_ex_cell{sub}.pre_cue_ex_cell & Dpost_cell_cats{sub}.allresp_cells];
    %use_cells = ~Dpost_ex_cell{sub}.pre_cue_ex_cell;
    PL_cue_resp_mag = [PL_cue_resp_mag, nanmean(Dpost_cell_resp{sub}.NR_Cue_resp_pos(:,[use_cells])) - nanmean(Dpost_cell_resp{sub}.NR_Cue_base_pos(:,[use_cells]))];
    PL_rew_resp_mag = [PL_rew_resp_mag, nanmean(Dpost_cell_resp{sub}.NR_Rew_resp_pos(:,[use_cells])) - nanmean(Dpost_cell_resp{sub}.NR_Rew_base_pos(:,[use_cells]))];
end

%% plotting
figure;
subplot(1,2,1);
scatter(D1_rew_resp_mag, D1_cue_resp_mag); hold on;
errorbarxy(mean(D1_rew_resp_mag), mean(D1_cue_resp_mag), std(D1_rew_resp_mag)/sqrt(length(D1_rew_resp_mag)), std(D1_cue_resp_mag)/sqrt(length(D1_cue_resp_mag)));
ylim([-0.2 1]); xlim([-0.2 1]);
hline(0, 'k'); vline(0, 'k');
xlabel('Reward response magniture (df/f)');
ylabel('Cue response magniture (df/f)');
title(['Day 1 n=', num2str(length(D1_cue_resp_mag))]);

subplot(1,2,2);
scatter(PL_rew_resp_mag, PL_cue_resp_mag);hold on;
errorbarxy(mean(PL_rew_resp_mag), mean(PL_cue_resp_mag), std(PL_rew_resp_mag)/sqrt(length(PL_rew_resp_mag)), std(PL_cue_resp_mag)/sqrt(length(PL_cue_resp_mag)));
ylim([-0.2 1]); xlim([-0.2 1]);
hline(0, 'k'); vline(0, 'k');
xlabel('Reward response magniture (df/f)');
ylabel('Cue response magniture (df/f)');
title(['Post Learning n=', num2str(length(PL_rew_resp_mag))]);
suptitle('Responsive Cells');





