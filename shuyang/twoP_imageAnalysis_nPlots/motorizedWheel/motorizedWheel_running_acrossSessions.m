%% assign document paths and experimental sessions
clear;
sessions = {'200116_img1041','200214_img1042','200217_img1061','200225_img1049',...
    '200319_img1064','200319_img1064_2'};
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\running\imaging_analysis\';
% behavior analysis results
days = {'1041-200116_1','1042-200114_1','1049-200225_1','1061-200217_1',...
    '1064-200319_1','1064-200319_2'};
color_code = {'b','r','k','c'};

%% df/f runtrig ave and runoffset ave and speed transit
dfOvF_runoff_fast_cells_allsession = []; %frame*cell
dfOvF_runoff_slow_cells_allsession = [];
dfOvF_runtrig_fast_cells_allsession = [];
dfOvF_runtrig_slow_cells_allsession = [];
dfOvF_spd_decrease_cells_allsession = [];
dfOvF_spd_increase_cells_allsession = [];

nPCs_total = 0;

for ii = 1:length(sessions)
    %load things
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvF_runoff_fast_cells = dfOvF_output.dfOvF_runoff_fast_cells; %frames*cells
    dfOvF_runoff_slow_cells = dfOvF_output.dfOvF_runoff_slow_cells; %frames*cells
    dfOvF_runtrig_fast_cells = dfOvF_output.dfOvF_runtrig_fast_cells;
    dfOvF_runtrig_slow_cells = dfOvF_output.dfOvF_runtrig_slow_cells;
    dfOvF_spd_decrease_cells = dfOvF_output.dfOvF_spd_decrease_cells;
    dfOvF_spd_increase_cells = dfOvF_output.dfOvF_spd_increase_cells;
    
    nPCs_total = nPCs_total + size(dfOvF_runoff_fast_cells,2);
    
    dfOvF_runoff_fast_cells_allsession = cat(2,dfOvF_runoff_fast_cells_allsession,dfOvF_runoff_fast_cells);%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
    dfOvF_runoff_slow_cells_allsession = cat(2,dfOvF_runoff_slow_cells_allsession,dfOvF_runoff_slow_cells);
    dfOvF_runtrig_fast_cells_allsession = cat(2,dfOvF_runtrig_fast_cells_allsession,dfOvF_runtrig_fast_cells);
    dfOvF_runtrig_slow_cells_allsession = cat(2,dfOvF_runtrig_slow_cells_allsession,dfOvF_runtrig_slow_cells);
    dfOvF_spd_decrease_cells_allsession = cat(2,dfOvF_spd_decrease_cells_allsession,dfOvF_spd_decrease_cells);
    dfOvF_spd_increase_cells_allsession = cat(2,dfOvF_spd_increase_cells_allsession,dfOvF_spd_increase_cells); 
end
        
x = (1:size(dfOvF_spd_increase_cells,1))/30;
dfOvF_across = figure;
subplot(2,3,1);
fast_errbar(x,dfOvF_runtrig_slow_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylabel('df/f'); ylim([0.05 1.05]); xlim([0 2.5]);
vline(1,'k');
title('stay-slow speed');
%xlim([0,1]);
subplot(2,3,4);
fast_errbar(x,dfOvF_runoff_slow_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylabel('df/f'); ylim([0.05 1.05]);xlim([0 2.5]);
vline(1,'k');
title('slow speed-stay');
subplot(2,3,2);
fast_errbar(x,dfOvF_runtrig_fast_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0.05 1.05]); vline(1,'k'); xlabel('time(s)');xlim([0 2.5]);
title('stay-fast speed');
subplot(2,3,5);
fast_errbar(x,dfOvF_runoff_fast_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0.05 1.05]); vline(1,'k'); xlabel('time(s)');xlim([0 2.5]);
title('fast speed-stay');
subplot(2,3,3);
fast_errbar(x,dfOvF_spd_decrease_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0.05 1.05]); xlim([0 2.5]);vline(1,'k');
title('speed decrease');
subplot(2,3,6);
fast_errbar(x,dfOvF_spd_increase_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
xlim([0 2.5]); ylim([0.05 1.05]); vline(1,'k');
title('speed increase');
supertitle(['dfOvF across sessions nPCs = ' num2str(nPCs_total)]);

savefig(['Z:\Analysis\motorizedWheel_Analysis\running\across_sessions\' 'across_sessions_dfOvF.fig']);

%% FR run trig ave and runoff ave

FR_runoff_fast_cells_allsession = []; %frame*cell
FR_runoff_slow_cells_allsession = [];
FR_runtrig_fast_cells_allsession = [];
FR_runtrig_slow_cells_allsession = [];
FR_spd_decrease_cells_allsession = [];
FR_spd_increase_cells_allsession = [];

nPCs_total = 0;

for ii = 1:length(sessions)
    %load things
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    spk_output = load([image_analysis_dest sessions{ii},'_spk_deconvolve_threshold-4.mat']);
    FR_runoff_fast_cells = spk_output.FR_runoff_fast_cells; %frames*cells
    FR_runoff_slow_cells = spk_output.FR_runoff_slow_cells; %frames*cells
    FR_runtrig_fast_cells = spk_output.FR_runtrig_fast_cells;
    FR_runtrig_slow_cells = spk_output.FR_runtrig_slow_cells;
    FR_spd_decrease_cells = spk_output.FR_spd_decrease_cells;
    FR_spd_increase_cells = spk_output.FR_spd_increase_cells;
    
    nPCs_total = nPCs_total + size(FR_runoff_fast_cells,2);
    
    FR_runoff_fast_cells_allsession = cat(2,FR_runoff_fast_cells_allsession,FR_runoff_fast_cells);%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
    FR_runoff_slow_cells_allsession = cat(2,FR_runoff_slow_cells_allsession,FR_runoff_slow_cells);
    FR_runtrig_fast_cells_allsession = cat(2,FR_runtrig_fast_cells_allsession,FR_runtrig_fast_cells);
    FR_runtrig_slow_cells_allsession = cat(2,FR_runtrig_slow_cells_allsession,FR_runtrig_slow_cells);
    FR_spd_decrease_cells_allsession = cat(2,FR_spd_decrease_cells_allsession,FR_spd_decrease_cells);
    FR_spd_increase_cells_allsession = cat(2,FR_spd_increase_cells_allsession,FR_spd_increase_cells); 
end
        
x = (1:size(FR_spd_increase_cells,1))/30;
FR_across = figure;
subplot(2,3,1);
fast_errbar(x,FR_runtrig_slow_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylabel('firing rate'); ylim([0 5.5]); xlim([0 2.5]);
vline(1,'k');
title('stay-slow speed');
%xlim([0,1]);
subplot(2,3,4);
fast_errbar(x,FR_runoff_slow_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylabel('firing rate'); ylim([0 5.5]);xlim([0 2.5]);
vline(1,'k');
title('slow speed-stay');
subplot(2,3,2);
fast_errbar(x,FR_runtrig_fast_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0 5.5]); vline(1,'k'); xlabel('time(s)');xlim([0 2.5]);
title('stay-fast speed');
subplot(2,3,5);
fast_errbar(x,FR_runoff_fast_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0 5.5]); vline(1,'k'); xlabel('time(s)');xlim([0 2.5]);
title('fast speed-stay');
subplot(2,3,3);
fast_errbar(x,FR_spd_decrease_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
ylim([0 5.5]); xlim([0 2.5]);vline(1,'k');
title('speed decrease');
subplot(2,3,6);
fast_errbar(x,FR_spd_increase_cells_allsession,2,'color',[0.1373 0.5451 0.2706],'shaded',true);
xlim([0 2.5]); ylim([0 5.5]); vline(1,'k');
title('speed increase');
supertitle(['FR across sessions nPC = ' num2str(nPCs_total)]);

savefig(['Z:\Analysis\motorizedWheel_Analysis\running\across_sessions\' 'across_sessions_FR.fig']);
