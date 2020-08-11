%This script is to

%% Section I: set paths and create analysis folders for each session
%define the directory and files
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'191114_img1040','191115_img1039','191115_img1041','191115_img1042','200316_img1064_airpuff_2'};
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
% behavior analysis results
days = {'1040-191114_1','1039-191115_1','1041-191115_1','1042-191115_1','1064-200316_2'};
color_code = {'c','r','y','g'};

%%  sectionII: airpuff triggered response-spikes
for ii = 1: length(sessions)
    % are there any spatial distribution based on response amplitude?
    % for each trial, the 16th frame is airpuff onset
    % consider 17-23th frame as the response window, max value is the response peak
    % calculate baseline F for each cell
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\']; 
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvFbtm_airpuff_stay_cells = dfOvF_strct.dfOvFbtm_airpuff_stay_cells;
    baseline_cell = zeros(1,size(dfOvFbtm_airpuff_stay_cells,2));
    resp_cell = zeros(1,size(dfOvFbtm_airpuff_stay_cells,2));
    for c = 1: size(dfOvFbtm_airpuff_stay_cells,2)
        % calculate baseline F for each cell
        temp = dfOvFbtm_airpuff_stay_cells(:,c);
        baseline_cell(c) = mean(temp(1:15));
        resp_cell(c) = max(temp(17:23)) - baseline_cell(c);
        %max_cell = max(temp);
        %max_inx = temp(temp == max_cell);
    end
    filename2 = dir([image_analysis_dest 'getTC\' '*' 'thresh97_coor0.8_mask3D.mat']);
    mask3D = load([image_analysis_dest 'getTC\' filename2.name]);
    mask3D = mask3D.mask3D;
    
    mask_cell = zeros(size(mask3D(:,:,1)));
    for c = 1: size(dfOvFbtm_airpuff_stay_cells,2)
        mask_cell = mask_cell + mask3D(:,:,c).*resp_cell(c);
    end
    figure; imagesc(mask_cell);
    title([sessions{ii} 'airpuff response amplitude of individual cells']);
    savefig([image_analysis_dest '\' sessions{ii} '_respAmp_mask']);
    
    % start with same y
    start_point = mean(dfOvFbtm_airpuff_stay_cells(1,:));
    
    figure;
    x1 = (1:size(dfOvFbtm_airpuff_stay_cells,1))/30;
    for c = 1: size(dfOvFbtm_airpuff_stay_cells,2)
        diff = dfOvFbtm_airpuff_stay_cells(1,c) - start_point;
        dfOvFbtm_airpuff_stay_cells_plot = dfOvFbtm_airpuff_stay_cells(:,c) - diff;
        plot(x1,dfOvFbtm_airpuff_stay_cells_plot); hold on;
    end
    xlim([0,1]);
    vline(16/30,'k');
    xlabel('time(s)');
    savefig([image_analysis_dest '\' sessions{ii} '_airresp_cells']);
end
    