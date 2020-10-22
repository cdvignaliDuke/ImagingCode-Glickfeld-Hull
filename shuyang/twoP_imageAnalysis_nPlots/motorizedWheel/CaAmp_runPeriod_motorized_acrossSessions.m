% across sessions: how is the Ca event amplitude different during the
% beginning of running, middle of running, and late running?
% also how many events are there during each running period

%% Section I: set paths and create analysis folders 
%define the directory and files
clear;
sessions = {'200116_img1041','200214_img1042','200217_img1061','200225_img1049','200319_img1064','200319_img1064_2'};
days = {'1041-200116_1','1042-200214_1','1061-200217_1','1049-200225_1','1064-200319_1','1064-200319_2'};
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\running\imaging_analysis\'; 
color_code = {'c','r','y','g'};

%% 

dfOvF_spk_run_begin_neurons_all = [];
dfOvF_spk_run_middle_neurons_all = [];
dfOvF_spk_run_late_neurons_all = [];
isospk_run_begin_all = 0;
isospk_run_middle_all = 0;
isospk_run_late_all = 0;
nPCs_Allsessions = zeros(1,length(sessions));

for ii = 1: length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    isoCa_output = load([image_analysis_dest sessions{ii} '_isoCaEvent.mat']);
    dfOvF_spk_run_begin_neurons_all = cat(1,dfOvF_spk_run_begin_neurons_all,isoCa_output.dfOvF_spk_run_begin_neurons);
    dfOvF_spk_run_middle_neurons_all = cat(1,dfOvF_spk_run_middle_neurons_all,isoCa_output.dfOvF_spk_run_middle_neurons);
    dfOvF_spk_run_late_neurons_all = cat(1,dfOvF_spk_run_late_neurons_all,isoCa_output.dfOvF_spk_run_late_neurons);
    isospk_run_begin_all = isospk_run_begin_all + isoCa_output.isospk_run_begin;
    isospk_run_middle_all = isospk_run_middle_all + isoCa_output.isospk_run_middle;
    isospk_run_late_all = isospk_run_late_all + isoCa_output.isospk_run_late;
    
    % total number of neurons in each session
    nPCs_Allsessions(ii) = size(isoCa_output.dfOvF_spk_run_begin_neurons,1);  
   
end
ntotalPCs = sum(nPCs_Allsessions);
x = (1:30)/30;

% when does isolated event happen during running?
total_event = isospk_run_begin_all+isospk_run_middle_all+isospk_run_late_all;
p_isospk_run_begin_all = isospk_run_begin_all*100/total_event;
p_isospk_run_middle_all = isospk_run_middle_all*100/total_event;
p_isospk_run_late_all = isospk_run_late_all*100/total_event;

figure; 
bar([p_isospk_run_begin_all,p_isospk_run_middle_all,p_isospk_run_late_all]);
set(gca,'XTickLabel',{'0-0.5s','0.5-1.5s','later than 1.5s'});
ylabel('percentage (%)');
title('time of isolated peaks during running windows');
savefig(['Z:\Analysis\motorizedWheel_Analysis\running\across_sessions\' 'across_sessions_isoCaEvent_run_time.fig']);

% how does time of event during running influecne amplitude
figure;
fast_errbar(x,dfOvF_spk_run_begin_neurons_all,1,'color',[0.1373 0.5451 0.2706]);hold on;
fast_errbar(x,dfOvF_spk_run_middle_neurons_all,1,'color',[0.2549 0.6706 0.3647]); hold on;
fast_errbar(x,dfOvF_spk_run_late_neurons_all,1,'color',[0.7294 0.8941 0.7020]); hold on;
legend('0-0.5s','0.5-1s','later than 1s');
text(0.75,0.3,['n total PCs = ' num2str(sum(ntotalPCs))]);
title('df/f of isolated events during different time period of running trials');
ylabel('df/f');
xlabel('time(s)');
savefig(['Z:\Analysis\motorizedWheel_Analysis\running\across_sessions\' 'across_sessions_CaAmpRun3conditions.fig']);

