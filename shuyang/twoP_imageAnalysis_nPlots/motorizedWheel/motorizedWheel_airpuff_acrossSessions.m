%% assign document paths and experimental sessions
clear;
sessions = {'200305_img1049','200319_img1064_airpuff','200319_img1064_airpuff_2'};
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\airpuff\imaging_analysis\';
% behavior analysis results
days = {'1049-200305_1','1064-200319_1','1064-200319_2'};
color_code = {'b','r','k','c'};

%% df/f 
dfOvF_airpuff_fast_cells_allsession = []; %frame*cell
dfOvF_airpuff_slow_cells_allsession = [];
dfOvF_airpuff_stay_cells_allsession = [];

nPCs_total = 0;

for ii = 1:length(sessions)
    %load things
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvF_airpuff_fast_cells = dfOvF_output.dfOvF_airpuff_fast_cells; %frames*cells
    dfOvF_airpuff_slow_cells = dfOvF_output.dfOvF_airpuff_slow_cells; %frames*cells
    dfOvF_airpuff_stay_cells = dfOvF_output.dfOvF_airpuff_stay_cells;
    
    nPCs_total = nPCs_total + size(dfOvF_airpuff_fast_cells,2);
    
    dfOvF_airpuff_fast_cells_allsession = cat(2,dfOvF_airpuff_fast_cells_allsession,dfOvF_airpuff_fast_cells);%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
    dfOvF_airpuff_slow_cells_allsession = cat(2,dfOvF_airpuff_slow_cells_allsession,dfOvF_airpuff_slow_cells);
    dfOvF_airpuff_stay_cells_allsession = cat(2,dfOvF_airpuff_stay_cells_allsession,dfOvF_airpuff_stay_cells);
end
        
x = (1:size(dfOvF_airpuff_fast_cells,1))/30;
dfOvF_across = figure;
subplot(1,3,1);
fast_errbar(x,dfOvF_airpuff_fast_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylabel('df/f'); ylim([0.05 1.05]); xlim([0 2]);
vline(1,'k');
title('fast speed');
%xlim([0,1]);
subplot(1,3,2);
fast_errbar(x,dfOvF_airpuff_slow_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylabel('df/f'); ylim([0.05 1.05]);xlim([0 2]);
vline(1,'k');
title('slow speed');
subplot(1,3,3);
fast_errbar(x,dfOvF_airpuff_stay_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0.05 1.05]); vline(1,'k'); xlabel('time(s)');xlim([0 2]);
title('staytionary');

supertitle(['tactile stimulus nPCs = ' num2str(nPCs_total)]);

savefig(['Z:\Analysis\motorizedWheel_Analysis\airpuff\across_sessions\' 'across_sessions_dfOvF_airpuff_motorized.fig']);

%% FR run trig ave and runoff ave

FR_airpuff_fast_cells_allsession = []; %frame*cell
FR_airpuff_slow_cells_allsession = [];
FR_airpuff_stay_cells_allsession = [];

nPCs_total = 0;

for ii = 1:length(sessions)
    %load things
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    spk_output = load([image_analysis_dest sessions{ii},'_spk_deconvolve_threshold-4.mat']);
    FR_airpuff_fast_cells = spk_output.FR_airpuff_fast_cells; %frames*cells
    FR_airpuff_slow_cells = spk_output.FR_airpuff_slow_cells; %frames*cells
    FR_airpuff_stay_cells = spk_output.FR_airpuff_stay_cells;
    
    nPCs_total = nPCs_total + size(FR_airpuff_fast_cells,2);
    
    FR_airpuff_fast_cells_allsession = cat(2,FR_airpuff_fast_cells_allsession,FR_airpuff_fast_cells);%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
    FR_airpuff_slow_cells_allsession = cat(2,FR_airpuff_slow_cells_allsession,FR_airpuff_slow_cells);
    FR_airpuff_stay_cells_allsession = cat(2,FR_airpuff_stay_cells_allsession,FR_airpuff_stay_cells);
end
        
x = (1:size(FR_airpuff_fast_cells,1))/30;
FR_across = figure;
subplot(1,3,1);
fast_errbar(x,FR_airpuff_fast_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylabel('df/f'); ylim([0.05 8]); xlim([0 2]);
vline(1,'k');
title('fast speed');
%xlim([0,1]);
subplot(1,3,2);
fast_errbar(x,FR_airpuff_slow_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylabel('df/f'); ylim([0.05 8]);xlim([0 2]);
vline(1,'k');
title('slow speed');
subplot(1,3,3);
fast_errbar(x,FR_airpuff_stay_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0.05 8]); vline(1,'k'); xlabel('time(s)');xlim([0 2]);
title('staytionary');

supertitle(['tactile stimulus nPCs = ' num2str(nPCs_total)]);

savefig(['Z:\Analysis\motorizedWheel_Analysis\airpuff\across_sessions\' 'across_sessions_FR_airpuff_motorized.fig']);
