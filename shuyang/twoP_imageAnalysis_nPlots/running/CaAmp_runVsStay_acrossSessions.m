% use this script after "CaAmp_runVsStayAirpuff"
% this script is comparing amplitude of well isolated Ca transients of all
% neurons in all sessions
% does the amplitude of Ca transient change under different conditions?
% condition 1 : stationary (not including 1s after running offset and 1s before running onset)
% condition 2: running (not including 1s after running onset and 1s? before running offset)
% condition 3: airpuff response during running (within 200ms of airpuff onset)
% condition 4: airpuff response during stationary (within 200ms of airpuff onset)
% events need to be isolated: no other events within 600ms before and/or after
% frm_stay_midpart: all stationary frames without beginning and end of stationary
% isolated_inx_neurons: index of isolated event peaks 1*#of cells, each cell element is a vector
% spk_iso_stay_neurons: index of isolated event peaks during stationary, 1*# of neurons, each cell element is a vector
% isoevent_stay_neurons: index of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
% dfOvF_btm_stay_isoevents: %df/f of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
% dfOvF_btm_stay_isoevents_neurons: average df/f of isolated event of each cell, averaged across events, neurons*frames
% ave_dfOvF_btm_stay_iso_session: average df/f of isolated event across cells, use this to plot

%% Section I: set paths and create analysis folders 
%define the directory and files
clear;
sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};
%sessions = {'190117_img1016'};
%days = {'1016-190117_1'};

image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; 
color_code = {'c','r','y','g'};

%% SectionII: for each session: running experiment, event amplitude during running vs stationary
% get events from each session, and save all of them into one variable at
% the end of the for loop
% for ii = 1: length(sessions)
%    % load data
%     image_analysis_dest = [image_analysis_base, sessions{ii}, '\']; 
%     dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
%     dfOvF_btm = dfOvF_strct.dfOvF_btm_cl;
%     behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{ii} '\'];
%     behav_output = load([behav_dest days{ii} '_behavAnalysis.mat']);
%     %behav_output = load([behav_dest days{i} '_first18000frames_behavAnalysis.mat']);
%     frm_stay_cell = behav_output.frames_stay_cell;
%     frm_run_cell = behav_output.frames_run_cell;
%     % stationary: not including 1s (30 frames) after running offset and 1s before
%     frm_stay_midpart_cell = {};
%     for t = 1: size(frm_stay_cell,2)
%         if length(frm_stay_cell{t})>= 61
%         frm_stay_midpart_cell = [frm_stay_midpart_cell frm_stay_cell{t}(31:end-30)];
%         end
%     end
%     frm_stay_midpart = cell2mat(frm_stay_midpart_cell); % all stationary frames without beginning and end of stationary
%     
%     %running (not including 1s after running onset and 1s? before running offset)
%     frm_run_midpart_cell = {};
%     for t = 1: size(frm_run_cell,2)
%         if length(frm_run_cell{t})>= 61
%         frm_run_midpart_cell = [frm_run_midpart_cell frm_run_cell{t}(31:end-30)];
%         end
%     end
%     if isempty(frm_run_midpart_cell) == 0
%         frm_run_midpart = cell2mat(frm_run_midpart_cell);
%     end
%     
%     % for each cell,find isolated events in running and stationary (no event within 650ms before that event)
%     threshold = -4;
%     spk_deconv_output = load([image_analysis_dest sessions{ii},'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
%     spk_inx_neurons = spk_deconv_output.spk_inx_cl;
%     
%     % identify events that are isolated: 600ms = 18 frames need to do this then see what behavioral state they're in, not the other way because
%     % if you find intersect of peak index and behavioral state frame first, then the diff between the frames you find can be big if there's a different behavioral state in between.
%     noother_after_cells = cell(1,size(spk_inx_neurons,2));
%     noother_befo_cells = cell(1,size(spk_inx_neurons,2));
%     isolated_inx_neurons = cell(1,size(spk_inx_neurons,2));
%     for c = 1: size(spk_inx_neurons,2)
%         noother_after_inx = find(diff(spk_inx_neurons{c})>18);% this means the indexed frame# being found here has a difference of >18 with the next frame#
%         noother_after_cells{c} = spk_inx_neurons{c}(noother_after_inx);
%         noother_befo_inx = noother_after_inx+1; % the diff of inx-(inx-1) > 18
%         noother_befo_cells{c} = spk_inx_neurons{c}(noother_befo_inx);
%         isolated_inx_neurons{c} = intersect(noother_after_cells{c},noother_befo_cells{c}); %index of isolated event peaks 1*#of cells, each cell element is a vector
%     end
%     
%     %isolated event peaks and matrices during stationary and running
%     spk_iso_stay_neurons = cell(1,size(spk_inx_neurons,2));% index of isolated event peaks during stationary, 1*# of neurons, each cell element is a vector
%     spk_iso_run_neurons = cell(1,size(spk_inx_neurons,2));
%     isoevent_stay_neurons = cell(1, size(spk_inx_neurons,2));% index of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
%     isoevent_run_neurons = cell(1, size(spk_inx_neurons,2));
%     dfOvF_btm_stay_isoevents = cell(1, size(spk_inx_neurons,2)); %df/f of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
%     dfOvF_btm_run_isoevents = cell(1, size(spk_inx_neurons,2));
%     dfOvF_btm_stay_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);% average df/f of isolated event of each cell, averaged across events, neurons*frames
%     dfOvF_btm_run_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);
%     
%     for c = 1: size(spk_inx_neurons,2)
%         %stationary
%         spk_iso_stay_neurons{c} = intersect(frm_stay_midpart,isolated_inx_neurons{c});
%         %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
%         isoevent_stay_neurons{c} = zeros(length(spk_iso_stay_neurons{c}),30); % for each cell, create a matrix, # of events*30(1s), events*frames, each row is a stationary event
%         dfOvF_btm_stay_isoevents{c} =  zeros(length(spk_iso_stay_neurons{c}),30);
%         dfOvF_c = dfOvF_btm(:,c);
%         for s = 1:length(spk_iso_stay_neurons{c})
%              isoevent_stay_neurons{c}(s,:) =  spk_iso_stay_neurons{c}(s)-11:spk_iso_stay_neurons{c}(s)+18;
%              dfOvF_btm_stay_isoevents{c}(s,:) = dfOvF_c(isoevent_stay_neurons{c}(s,:));
%         end
%         % average across trials
%         dfOvF_btm_stay_isoevents_neurons(c,:) = mean(dfOvF_btm_stay_isoevents{c}); 
%         %running
%         if isempty(frm_run_midpart_cell) == 0 % if there are running trials that are longer than 2s, then continue
%             spk_iso_run_neurons{c} = intersect(frm_run_midpart,isolated_inx_neurons{c});
%             %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
%             isoevent_run_neurons{c} = zeros(length(spk_iso_run_neurons{c}),30);
%             dfOvF_btm_run_isoevents{c} =  zeros(length(spk_iso_run_neurons{c}),30);
%             for s = 1:length(spk_iso_run_neurons{c})
%                 isoevent_run_neurons{c}(s,:) =  spk_iso_run_neurons{c}(s)-11:spk_iso_run_neurons{c}(s)+18;
%                 dfOvF_btm_run_isoevents{c}(s,:) = dfOvF_c(isoevent_run_neurons{c}(s,:));
%             end
%             % average across trials
%             dfOvF_btm_run_isoevents_neurons(c,:) = mean(dfOvF_btm_run_isoevents{c});
%         end
%     end
%     
%     % average across cells
%     
%     ave_dfOvF_btm_stay_iso_session = mean(dfOvF_btm_stay_isoevents_neurons);
%     % ste
%     ste_dfOvF_btm_stay_iso_session = std(dfOvF_btm_stay_isoevents_neurons,0,1)/sqrt(size(dfOvF_btm_stay_isoevents_neurons,1));
%     
%     if isempty(frm_run_midpart_cell) == 0
%         % if there's a cell that doesn't spike during running midpart,there will be a line of NaNs in dfOvF_btm_run_isoevents_neurons, need to delete this Nan line before doing average
%         dfOvF_btm_run_isoevents_neurons = dfOvF_btm_run_isoevents_neurons(all(~isnan(dfOvF_btm_run_isoevents_neurons),2),:);
%         ave_dfOvF_btm_run_iso_session = mean(dfOvF_btm_run_isoevents_neurons);
%         ste_dfOvF_btm_run_iso_session = std(dfOvF_btm_run_isoevents_neurons,0,1)/sqrt(size(dfOvF_btm_stay_isoevents_neurons,1));
%     end
%     
%     %count how many events you have in each behavioral condition and how many cells in each session 
%     nevents_run = sum(cellfun(@(x)x(1),cellfun(@size,isoevent_run_neurons,'UniformOutput',false)));% # of isolated events for all cells during running
%     nevents_stay = sum(cellfun(@(x)x(1),cellfun(@size,isoevent_stay_neurons,'UniformOutput',false)));
%     ncells = size(spk_inx_neurons,2);
%     %plot
%     x = (1:30)/30;
%     Caevent = figure;
%     errorbar(x,ave_dfOvF_btm_stay_iso_session,ste_dfOvF_btm_stay_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
%     errorbar(x,ave_dfOvF_btm_run_iso_session,ste_dfOvF_btm_run_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
%     % number of calcium events is wrong in the current legend
%     %legend(['stationary nevents=' num2str(nevents_stay)],['running nevents=' num2str(nevents_run)]);
%     legend(['nPCs= ' num2str(ncells)]);
%     ylabel('df/f');
%     xlabel('time(s)');
%     title(sessions(ii));
%     
%     saveas(Caevent,[image_analysis_dest '\' days{ii} '_CaAmp']);
%     % save variables
%     
%     if isempty(frm_run_midpart_cell) == 0
%         save([image_analysis_dest sessions{ii} '_isoCaEvent.mat'],'frm_stay_midpart',...
%             'frm_run_midpart','isolated_inx_neurons','spk_iso_stay_neurons','spk_iso_run_neurons',...
%             'isoevent_stay_neurons','isoevent_run_neurons','dfOvF_btm_stay_isoevents','dfOvF_btm_run_isoevents',...
%             'dfOvF_btm_stay_isoevents_neurons','dfOvF_btm_run_isoevents_neurons',...
%             'ave_dfOvF_btm_stay_iso_session','ste_dfOvF_btm_stay_iso_session',...
%             'ave_dfOvF_btm_run_iso_session','ste_dfOvF_btm_run_iso_session','nevents_run','nevents_stay','ncells');
%     else
%         fprintf(['no running trials longer than 2s in this session ' sessions{ii}]);
%         save([image_analysis_dest sessions{ii} '_isoCaEvent.mat'],'frm_stay_midpart',...
%             'frm_run_midpart','isolated_inx_neurons','spk_iso_stay_neurons',...
%             'isoevent_stay_neurons','dfOvF_btm_stay_isoevents','dfOvF_btm_stay_isoevents_neurons',...
%             'ave_dfOvF_btm_stay_iso_session','ste_dfOvF_btm_stay_iso_session','nevents_stay');
%     end
% end  

%% 
ave_dfOvF_btm_stay_iso_Allsessions = zeros(length(sessions),30);
ste_dfOvF_btm_stay_iso_Allsessions = zeros(length(sessions),30);
ave_dfOvF_btm_run_iso_Allsessions = zeros(length(sessions),30);
ste_dfOvF_btm_run_iso_Allsessions = zeros(length(sessions),30);
nevents_stay_Allsessions = zeros(1,length(sessions));
nevents_run_Allsessions = zeros(1,length(sessions));
nPCs_Allsessions = zeros(1,length(sessions));
nPCs_run = zeros(1,length(sessions));%cells included when comparing amplitude during running. only cells fire during early, middle, abd late running are included because we need to draw the same number of events from every condition
dfOvF_btm_run_isoevents_neurons_all = [];
dfOvF_btm_stay_isoevents_neurons_all = [];
dfOvF_spk_run_begin_neurons_all = [];
dfOvF_spk_run_middle_neurons_all = [];
dfOvF_spk_run_late_neurons_all = [];
isospk_run_begin_all = 0;
isospk_run_middle_all = 0;
isospk_run_late_all = 0;
for ii = 1: length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    isoCa_output = load([image_analysis_dest sessions{ii} '_isoCaEvent.mat']);
    ave_dfOvF_btm_stay_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_stay_iso_session;
    ave_dfOvF_btm_run_iso_Allsessions(ii,:) = isoCa_output.ave_dfOvF_btm_run_iso_session;
    dfOvF_btm_run_isoevents_neurons_all = cat(1,dfOvF_btm_run_isoevents_neurons_all,isoCa_output.dfOvF_btm_run_isoevents_neurons);
    dfOvF_btm_stay_isoevents_neurons_all = cat(1,dfOvF_btm_stay_isoevents_neurons_all,isoCa_output.dfOvF_btm_stay_isoevents_neurons);
    dfOvF_spk_run_begin_neurons_all = cat(1,dfOvF_spk_run_begin_neurons_all,isoCa_output.dfOvF_spk_run_begin_neurons);
    dfOvF_spk_run_middle_neurons_all = cat(1,dfOvF_spk_run_middle_neurons_all,isoCa_output.dfOvF_spk_run_middle_neurons);
    dfOvF_spk_run_late_neurons_all = cat(1,dfOvF_spk_run_late_neurons_all,isoCa_output.dfOvF_spk_run_late_neurons);
    nPCs_run (ii) = size(isoCa_output.dfOvF_spk_run_begin_neurons,1);
    isospk_run_begin_all = isospk_run_begin_all + isoCa_output.isospk_run_begin;
    isospk_run_middle_all = isospk_run_middle_all + isoCa_output.isospk_run_middle;
    isospk_run_late_all = isospk_run_late_all + isoCa_output.isospk_run_late;

    % total number of events from all neurons in each session
%     nevents_stay_Allsessions(ii) = isoCa_output.nevents_stay;
%     nevents_run_Allsessions(ii) = isoCa_output.nevents_run;
    
    % total number of neurons in each session
    nPCs_Allsessions(ii) = size(isoCa_output.dfOvF_btm_stay_isoevents_neurons,1);  
   
end
% average across sessions
ave_dfOvF_btm_stay_iso_across = mean(ave_dfOvF_btm_stay_iso_Allsessions);
ave_dfOvF_btm_run_iso_across = mean(ave_dfOvF_btm_run_iso_Allsessions);

% need to change this later
ste_dfOvF_btm_stay_iso_across = std(dfOvF_btm_stay_isoevents_neurons_all)/sqrt(size(dfOvF_btm_stay_isoevents_neurons_all,1));
ste_dfOvF_btm_run_iso_across = std(dfOvF_btm_run_isoevents_neurons_all)/sqrt(size(dfOvF_btm_run_isoevents_neurons_all,1));

% ntotalevents_stay = sum(nevents_stay_Allsessions);
% ntotalevents_run = sum(nevents_run_Allsessions);
ntotalPCs = sum(nPCs_Allsessions);
x = (1:30)/30;
figure;
%haven't figured out how to use RGB colors in errorbar
errorbar(x,ave_dfOvF_btm_stay_iso_across,ste_dfOvF_btm_stay_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
errorbar(x,ave_dfOvF_btm_run_iso_across,ste_dfOvF_btm_run_iso_across,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
%legend(['stationary nevents=' num2str(ntotalevents_stay)],['running nevents=' num2str(ntotalevents_run)]);
legend('stataionry','running');
text(0.75,0.45, ['n total PCs=' num2str(ntotalPCs)]);
ylim([0 0.6]);
xlabel('time(s)'); ylabel('df/f');
savefig(['Z:\Analysis\2P_MovingDots_Analysis\across_sessions\' 'across_sessions_CaAmpRunvsStay.fig']);

% when does isolated event happen during running?
total_event = isospk_run_begin_all+isospk_run_middle_all+isospk_run_late_all;
p_isospk_run_begin_all = isospk_run_begin_all*100/total_event;
p_isospk_run_middle_all = isospk_run_middle_all*100/total_event;
p_isospk_run_late_all = isospk_run_late_all*100/total_event;

figure; 
bar([p_isospk_run_begin_all,p_isospk_run_middle_all,p_isospk_run_late_all]);
set(gca,'XTickLabel',{'0-0.5s','0.5-1s','later than 1s'});
ylabel('percentage (%)');
title('time of isolated peaks during running windows');
savefig(['Z:\Analysis\2P_MovingDots_Analysis\across_sessions\' 'across_sessions_isoCaEvent_run_time.fig']);

% how does time of event during running influecne amplitude
figure;
fast_errbar(x,dfOvF_spk_run_begin_neurons_all,1,'color',[0.1373 0.5451 0.2706]);hold on;
fast_errbar(x,dfOvF_spk_run_middle_neurons_all,1,'color',[0.2549 0.6706 0.3647]); hold on;
fast_errbar(x,dfOvF_spk_run_late_neurons_all,1,'color',[0.7294 0.8941 0.7020]); hold on;
legend('0-0.5s','0.5-1s','later than 1s');
text(0.75,0.3,['n total PCs = ' num2str(sum(nPCs_run))]);
title('df/f of isolated events during different time period of running trials');
ylabel('df/f');
xlabel('time(s)');
savefig(['Z:\Analysis\2P_MovingDots_Analysis\across_sessions\' 'across_sessions_CaAmpRun3conditions.fig']);

