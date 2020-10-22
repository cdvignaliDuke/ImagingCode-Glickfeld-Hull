%% 2 photon: dfOvF
% have to plot these after deconvolution. Because deconvolution tells you
% which cells are bad, and you don't look at the fluorescence of those cells (in section one, it's TCave_cl instead of TCave).
%plot average dfOvF for all cells, with speed.
%can comment out the speed part if don't want to plot speed

%% assign document paths and experimental sessions
clear;
sessions = {'191220_img1042'}; 
days = {'1042-191220_1'};
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
%image_analysis_base    = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
color_code = {'b','r','k','c'};

%% SECTION I df/f
for i = 1:size(sessions,2)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{i} '\'];
    %behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{i} '\'];  
    threshold = -4;
    filename = dir([image_analysis_dest '*' num2str(threshold) '_TCave_cl.mat']);
    %filename = dir([image_analysis_dest '*' 'first18000frames' '*' num2str(threshold) '_TCave_cl.mat']);
    TCave_cl = load([image_analysis_dest filename.name]);
    TCave_cl = TCave_cl.TCave_cl;
    behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
    %behav_output = load([behav_dest days{i} '_first18000frames_behavAnalysis.mat']);
    speed = double(behav_output.speed);
    frm_stay_cell = behav_output.frames_stay_cell;
    frm_stay = cell2mat(frm_stay_cell);

    %calculate baseline(F) and df/F
    % two ways to calculate baselineF: 
    % 1.stationary without airpuff (get rid of frames 300ms after airpuff onset)
    baseline = mean(TCave_cl(frm_stay,:));         % baseline for each cell
    dfOvF_stay = zeros(size(TCave_cl,1),size(TCave_cl,2));
    for n = 1:size(TCave_cl,2)                     % for each cell
        dfOvF_stay(:,n) = (TCave_cl(:,n)-baseline(n))/baseline(n); %frames*cell
    end
    
     % 2.bottom 10% of all fluorescence
    dfOvF_btm = zeros(size(TCave_cl,1),size(TCave_cl,2));
    baseline_btm = zeros(1, size(TCave_cl,2));
    for n = 1:size(TCave_cl,2)                     % for each cell
        TCave_sort = sort(TCave_cl(:,n),1,'ascend');
        btm = TCave_sort(1:floor(length(TCave_sort)*0.1));
        %        % test if btm always belongs to running
        %         btm_inx = frames(ismember(TCave_cl(:,n),btm));
        %         NbtmInRun = sum((ismember(frm_run,btm_inx))==1);
        baseline_btm(n) = mean(btm);
        dfOvF_btm(:,n) = (TCave_cl(:,n)-baseline_btm(n))/baseline_btm(n); %frames*cell
    end
    save([image_analysis_dest sessions{i} '_dfOvF.mat'],'dfOvF_stay','dfOvF_btm');
 
end