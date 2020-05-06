% use this script after "CaAmp_runVsStayAirpuff"
% this script is comparing amplitude of well isolated Ca transients of all
% neurons in all sessions
% does the amplitude of Ca transient change under different conditions?
% condition 1 : stationary (not including 1s after running offset and 1s before running onset)
% condition 2: fast running (not including 1s after running onset and 1s? before running offset)
% condition 3: slow running

%% Section I: set paths and create analysis folders 
%define the directory and files

clear;
sessions = {'200305_img1049','200319_img1064_airpuff','200319_img1064_airpuff_2'};
days = {'1049-200305_1','1064-200319_1','1064-200319_2'};
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\airpuff\imaging_analysis\';
color_code = {'c','r','y','g'};

%% 
ave_dfOvF_btm_stay_iso_Allsessions = zeros(length(sessions),30);
ave_dfOvF_btm_run_iso_Allsessions = zeros(length(sessions),30);
ave_dfOvF_btm_staypuff_iso_Allsessions = zeros(length(sessions),30);
nPCs_Allsessions = zeros(1,length(sessions));
dfOvF_btm_run_isoevents_neurons_all = [];
dfOvF_btm_stay_isoevents_neurons_all = [];
dfOvF_btm_staypuff_isoevents_neurons_all = [];

for ii = 1: length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    isoCa_output = load([image_analysis_dest sessions{ii} '_isoCaEvent.mat']);
    ave_dfOvF_btm_stay_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_stay_iso_session;
    ave_dfOvF_btm_run_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_run_iso_session;
    ave_dfOvF_btm_staypuff_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_staypuff_iso_session;
    dfOvF_btm_run_isoevents_neurons_all = cat(1,dfOvF_btm_run_isoevents_neurons_all,isoCa_output.dfOvF_btm_run_isoevents_neurons);
    dfOvF_btm_stay_isoevents_neurons_all = cat(1,dfOvF_btm_stay_isoevents_neurons_all,isoCa_output.dfOvF_btm_stay_isoevents_neurons);
    dfOvF_btm_staypuff_isoevents_neurons_all = cat(1,dfOvF_btm_staypuff_isoevents_neurons_all, isoCa_output.dfOvF_btm_staypuff_isoevents_neurons);
    % total number of events from all neurons in each session
%     nevents_stay_Allsessions(ii) = isoCa_output.nevents_stay;
%     nevents_run_Allsessions(ii) = isoCa_output.nevents_run;
    
    % total number of neurons in each session
    nPCs_Allsessions(ii) = size(isoCa_output.dfOvF_btm_stay_isoevents_neurons,1);  
   
end
% average across sessions
ave_dfOvF_btm_stay_iso_across = mean(ave_dfOvF_btm_stay_iso_Allsessions);
ave_dfOvF_btm_run_iso_across = mean(ave_dfOvF_btm_run_iso_Allsessions);
ave_dfOvF_btm_staypuff_iso_across = mean(ave_dfOvF_btm_staypuff_iso_Allsessions);

ste_dfOvF_btm_stay_iso_across = std(dfOvF_btm_stay_isoevents_neurons_all)/sqrt(size(dfOvF_btm_stay_isoevents_neurons_all,1));
ste_dfOvF_btm_run_iso_across = std(dfOvF_btm_run_isoevents_neurons_all)/sqrt(size(dfOvF_btm_run_isoevents_neurons_all,1));
ste_dfOvF_btm_staypuff_iso_across = std(dfOvF_btm_staypuff_isoevents_neurons_all)/sqrt(size(dfOvF_btm_staypuff_isoevents_neurons_all,1));
% ntotalevents_stay = sum(nevents_stay_Allsessions);
% ntotalevents_run = sum(nevents_run_Allsessions);
ntotalPCs = sum(nPCs_Allsessions);
x = (1:30)/30;
figure;
%haven't figured out how to use RGB colors in errorbar
errorbar(x,ave_dfOvF_btm_stay_iso_across,ste_dfOvF_btm_stay_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
errorbar(x,ave_dfOvF_btm_run_iso_across,ste_dfOvF_btm_run_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
errorbar(x,ave_dfOvF_btm_staypuff_iso_across,ste_dfOvF_btm_staypuff_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
%legend(['stationary nevents=' num2str(ntotalevents_stay)],['running nevents=' num2str(ntotalevents_run)]);
legend('stationary', 'running','airpuff response during stationary');
text(0.75,0.7, ['n total PCs=' num2str(ntotalPCs)]);
%ylim([0 0.6]);
xlabel('time(s)'); ylabel('df/f');
savefig(['Z:\Analysis\motorizedWheel_Analysis\running\across_sessions' 'across_sessions_CaAmpRunvsStayVsAirpuff.fig']);

