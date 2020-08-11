%% 2 photon: dfOvF
% have to run this script after deconvolution. Because deconvolution tells you
% which cells are bad, and you don't look at the fluorescence of those cells (in section one, it's TCave_cl instead of TCave).
% using bottom 10% of fluorescence values as baseline

%% assign document paths and experimental sessions
clear;
sessions = {'200116_img1041','200214_img1042','200217_img1061','200225_img1049'}; 
% days = {'1042-191115'};
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\running\imaging_analysis\'; 
%image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\airpuff\imaging_analysis\'; 
%image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\tone\imaging_analysis\'; 
%image_analysis_base  = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
%image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder

color_code = {'b','r','k','c'};

%% SECTION I df/f
for i = 1:size(sessions,2)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    %behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\behavioral_analysis\' days{i} '\'];
    threshold = -4;
    filename = dir([image_analysis_dest '*' num2str(threshold) '_TCave_cl.mat']);
    TCave_cl = load([image_analysis_dest filename.name]);
    TCave_cl = TCave_cl.TCave_cl;
    
%     %if use all TCs
%     filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
%     TCave = load([image_analysis_dest 'getTC\' filename.name]);
%     TCave_cl = TCave.tc_avg;
%     
    %behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
    
    %calculate baseline(F) and df/F
    % bottom 10% of all fluorescence
    dfOvF_btm_cl = zeros(size(TCave_cl,1),size(TCave_cl,2));
    baseline_btm = zeros(1, size(TCave_cl,2));
    for n = 1:size(TCave_cl,2)                     % for each cell
        TCave_sort = sort(TCave_cl(:,n),1,'ascend');
        btm = TCave_sort(1:floor(length(TCave_sort)*0.1));
        %        % test if btm always belongs to running
        %         btm_inx = frames(ismember(TCave_cl(:,n),btm));
        %         NbtmInRun = sum((ismember(frm_run,btm_inx))==1);
        baseline_btm(n) = mean(btm);
        dfOvF_btm_cl(:,n) = (TCave_cl(:,n)-baseline_btm(n))/baseline_btm(n); %frames*cell
    end
    save([image_analysis_dest sessions{i} '_dfOvF.mat'],'dfOvF_btm_cl');
    
end