% Scropt for determining the fraction of responsive cells for various experimental conditions
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
    D1_ex_cell{sub} = load([dest_sub, '_pre_cue_ex_cell.mat']);
end

%% Post learning
Dpost_cell_cats = [];  Dpost_pvals = [];  Dpost_cell_resp = []; Dpost_ex_cell=[];
for sub = 1:size(days_post,2)
    %collect session/mouse information  set pathnames
    [dest_sub, ~, ~, ~, session, session_date, mouse_num, mouse_ID, rID] = get_sess_info(days_post{sub}, runID, crp_dir, []);

    %load data into cell
    Dpost_cell_cats{sub} = load([dest_sub '_cell_categories.mat']);
    Dpost_ex_cell{sub} = load([dest_sub, '_pre_cue_ex_cell.mat']);
end

%% Unexpected Reward
UR_cell_cats = [];  UR_pvals = [];  UR_cell_resp = []; UR_ex_cell=[];
for sub = 1:size(days_UR,2)
    %collect session/mouse information  set pathnames
    [dest_sub, ~, ~, ~, session, session_date, mouse_num, mouse_ID, rID] = get_sess_info(days_UR{sub}, runID, crp_dir, []);

    %load data into cell
    UR_cell_cats{sub} = load([dest_sub '_cell_categories.mat']);
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
    D1000_ex_cell{sub} = load([dest_sub, '_pre_cue_ex_cell.mat']);
end

%% collect proper data

%get number of cue vs rew responsive cells
D1_rew_frac =[];
PL_rew_frac =[];
D1_cue_frac = [];
PL_cue_frac = [];
for sub = 1:length(D1_cell_cats)
    if isempty(D1_cell_cats{sub}) 
        D1_rew_frac = [D1_rew_frac, NaN];
        D1_cue_frac = [D1_cue_frac, NaN];
    else
        %if sum([D1_cell_cats{sub}.allresp_cells & ~D1_ex_cell{sub}.pre_cue_ex_cell]) < 10  %if there are fewer than 10 resp cells
        if sum([~D1_ex_cell{sub}.pre_cue_ex_cell]) < 10  %if there are fewer than 10 cells
            D1_rew_frac = [D1_rew_frac, NaN];
            D1_cue_frac = [D1_cue_frac, NaN];
        else
            NR_rew_mat = zeros(1,D1_cell_cats{sub}.nCells); 
            NR_rew_mat([D1_cell_cats{sub}.NR_Rew_resp_cells_pos]) = 1;;
            NR_cue_mat = zeros(1,D1_cell_cats{sub}.nCells); 
            NR_cue_mat([D1_cell_cats{sub}.NR_Cue_resp_cells_pos]) = 1;;
%             D1_rew_frac = [D1_rew_frac, sum([NR_rew_mat & ~D1_ex_cell{sub}.pre_cue_ex_cell] ) / sum([D1_cell_cats{sub}.allresp_cells & ~D1_ex_cell{sub}.pre_cue_ex_cell])];
%             D1_cue_frac = [D1_cue_frac, sum([NR_cue_mat & ~D1_ex_cell{sub}.pre_cue_ex_cell ]) / sum([D1_cell_cats{sub}.allresp_cells & ~D1_ex_cell{sub}.pre_cue_ex_cell])];
            D1_rew_frac = [D1_rew_frac, sum([NR_rew_mat & ~D1_ex_cell{sub}.pre_cue_ex_cell] ) / sum([~D1_ex_cell{sub}.pre_cue_ex_cell])];
            D1_cue_frac = [D1_cue_frac, sum([NR_cue_mat & ~D1_ex_cell{sub}.pre_cue_ex_cell ]) / sum([~D1_ex_cell{sub}.pre_cue_ex_cell])];
            %===========================================
            
        end
    end
    if sum([Dpost_cell_cats{sub}.allresp_cells & ~Dpost_ex_cell{sub}.pre_cue_ex_cell]) < 10
        PL_rew_frac = [PL_rew_frac, NaN];
        PL_cue_frac = [PL_cue_frac, NaN];
    else
        NR_rew_mat = zeros(1,Dpost_cell_cats{sub}.nCells);
        NR_rew_mat([Dpost_cell_cats{sub}.NR_Rew_resp_cells_pos]) = 1;;
        NR_cue_mat = zeros(1,Dpost_cell_cats{sub}.nCells);
        NR_cue_mat([Dpost_cell_cats{sub}.NR_Cue_resp_cells_pos]) = 1;;
        PL_rew_frac = [PL_rew_frac, sum( [NR_rew_mat & ~Dpost_ex_cell{sub}.pre_cue_ex_cell] ) / sum([~Dpost_ex_cell{sub}.pre_cue_ex_cell])];
        PL_cue_frac = [PL_cue_frac, sum( [NR_cue_mat & ~Dpost_ex_cell{sub}.pre_cue_ex_cell] ) / sum([~Dpost_ex_cell{sub}.pre_cue_ex_cell])];
        %==============================================
        
    end
end

%find 1000ms overlap between cue and resp cells 
%find rew cells which are also cue cells    Also get num rew cells 
num_rew_cells = [];
num_cue_cells = [];
rew_and_cue_cells = [];
for sub = 1:length(D1000_cell_cats)
    if isempty(D1000_cell_cats{sub}) || sum(D1000_cell_cats{sub}.allresp_cells) < 10 
        num_rew_cells = [num_rew_cells, NaN];
        num_cue_cells = [num_cue_cells, NaN];
        rew_and_cue_cells = [rew_and_cue_cells, NaN];
    else
        num_rew_cells = [num_rew_cells, sum(D1000_cell_cats{sub}.rew_cells)];
        num_cue_cells = [num_cue_cells, sum(D1000_cell_cats{sub}.cue_cells)];
        rew_and_cue_cells = [rew_and_cue_cells, sum([D1000_cell_cats{sub}.rew_cells & D1000_cell_cats{sub}.cue_cells])];
    end
end

%% plotting 

%Day1 rew vs cue resp
figure;
subplot(1,2,1);
scatter(D1_rew_frac(isfinite(D1_rew_frac)), D1_cue_frac(isfinite(D1_cue_frac)), 'k'); hold on;
errorbarxy(mean(D1_rew_frac(isfinite(D1_rew_frac))), mean(D1_cue_frac(isfinite(D1_cue_frac))), std(D1_rew_frac(isfinite(D1_rew_frac)))/sqrt(sum(isfinite(D1_rew_frac))), std(D1_cue_frac(isfinite(D1_cue_frac)))/sqrt(sum(isfinite(D1_cue_frac))) );
hold on;
plot([0, 1], [0,1], 'k');
ylabel('fraction of cue resp neurons');
xlabel('fraction of rew resp neuron');
title('Day 1');
ylim([0 1.05]); xlim([0 1.05]);

%post learning
subplot(1,2,2);
scatter(PL_rew_frac(isfinite(PL_rew_frac)), PL_cue_frac(isfinite(PL_cue_frac)), 'k'); hold on;
errorbarxy(mean(PL_rew_frac(isfinite(PL_rew_frac))), mean(PL_cue_frac(isfinite(PL_cue_frac))), std(PL_rew_frac(isfinite(PL_rew_frac)))/sqrt(sum(isfinite(PL_rew_frac))), std(PL_cue_frac(isfinite(PL_cue_frac)))/sqrt(sum(isfinite(PL_cue_frac))) );
hold on;
plot([0, 1], [0,1], 'k');
ylabel('fraction of cue resp neurons');
xlabel('fraction of rew resp neuron');
title('Post Learning');
ylim([0 1.05]); xlim([0 1.05]);

% 1000ms delay % of cue resp also rew resp   vice versa
figure; 
% Fraction of CUE resp cells which are also REW resp
subplot(1,2,1);
these_fracs = [rew_and_cue_cells(isfinite(rew_and_cue_cells))./num_cue_cells(isfinite(num_cue_cells))];
plot(ones(1,length(these_fracs)), these_fracs, 'ok', 'LineStyle', 'none'); hold on;
errorbar(1, mean(these_fracs), std(these_fracs)/sqrt(length(these_fracs)), 'or');
xlabel('cue cells');
ylabel('fraction also reward responsive');
ylim([0 1]);  xlim([0 2]);
% Fraction of REW resp cells which are also CUE resp
subplot(1,2,2);
these_fracs = [rew_and_cue_cells(isfinite(rew_and_cue_cells))./num_rew_cells(isfinite(num_rew_cells))];
plot(ones(1,length(these_fracs)), these_fracs, 'ok', 'LineStyle', 'none'); hold on;
errorbar(1, mean(these_fracs), std(these_fracs)/sqrt(length(these_fracs)), 'or');
xlabel('Reward cells');
ylabel('fraction also cue responsive');
ylim([0 1]);    xlim([0 2]);
suptitle('1.1s delay condtion: fractional overlap between reward and cue responsive cells');

% figure;
% for sub = 1:length(D1_rew_frac)
%     plot([1,2],[D1_rew_frac(sub),PL_cue_frac(sub)], 'k'); hold on;
% end
% errorbar([1,2],[mean(D1_rew_frac),mean(PL_cue_frac)], [std(D1_rew_frac)/sqrt(length(D1_rew_frac)), std(PL_cue_frac)/sqrt(length(PL_cue_frac))], 'r', 'LineStyle', 'none');
% ylabel('fraction of responsive neurons');
% xlabel('Day 1 Reward resp                  post-learning cue resp');
% title('All responsive neurons (positive and negative)');
% ylim([0 1]);
% 
% figure; 
% for sub = 1:length(D1_cue_frac)
%     plot([1,2],[D1_cue_frac(sub),PL_cue_frac(sub)], 'k'); hold on;
% end
% errorbar([1,2],[mean(D1_cue_frac),mean(PL_cue_frac)], [std(D1_cue_frac)/sqrt(length(D1_cue_frac)), std(PL_cue_frac)/sqrt(length(PL_cue_frac))], 'r', 'LineStyle', 'none');
% ylabel('fraction of responsive neurons');
% xlabel('Day 1 cue resp                  post-learning cue resp');
% title('All responsive neurons (positive and negative)');
% ylim([0 1]);
% 
% figure;
% for sub = 1:length(D1_rew_frac)
%     plot([1,2],[D1_rew_frac(sub),PL_rew_frac(sub)], 'k'); hold on;
% end
% errorbar([1,2],[mean(D1_rew_frac),mean(PL_rew_frac)], [std(D1_rew_frac)/sqrt(length(D1_rew_frac)), std(PL_rew_frac)/sqrt(length(PL_rew_frac))], 'r', 'LineStyle', 'none');
% ylabel('fraction of responsive neurons');
% xlabel('Day 1 Reward resp                  post-learning rew resp');
% title('fraction of all neurons');
% ylim([0 1]);






