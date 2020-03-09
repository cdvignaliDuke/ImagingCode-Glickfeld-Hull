%This script is to plot the speed before and after airpuff onset. Speeds
%are aligned to stimulus onset. No analysis for imaging data
% make sure the criteria is what you want! The criteria is different from the
% script called airpuffStim_stateConsistent
%% Section I: set paths and create analysis folders for each session
%define the directory and files
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'190903_img1023_1'};
days = {'1023-190903_1'};
sessionID = {'1023-190903_1'};

%image_dest_base    = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'c','r','y','g'};
%%
for ii = 1: length(sessions)
    behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    speed = behav_struct.speed; speed = double(speed);
    frames_run_cell = behav_struct.frames_run_cell;
    frames_stay_cell = behav_struct.frames_stay_cell;
    airpuffon = behav_struct.airpuffon1;
    %image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    %dfOvF_strct = load([image_dest, '_dfOvF_staybase.mat']);
    %dfOvF = dfOvF_strct.dfOvF_staybase;
    
    before = 10;%# of frames
    after = 15;
    period = 30; %#of frames before and after airpuff onset
    
    % generate frame matrixes:
    % frms_airpuff_all: 30 frames (1or2s) before and 30 frames after stim, all stim trials together
    % frms_airtrig_stay: 10 or 20 frames (0.6s)stationary before stim and
    % 1s after stim (mice just needs to be )
    % frms_airtrig_run: 10 or 20 frames (0.6s) running before stim and 1s after stim
    
    frms_airpuff_all = zeros(period*2,length(airpuffon));
    
    for i = 1: length(airpuffon)
        if airpuffon(i)- period <1 || airpuffon(i)+period-1 > length(speed)
            continue
        else
            frms_airpuff_all(:,i) = airpuffon(i)-period:airpuffon(i)+period-1; % frames*trial
        end
    end
    frms_airpuff_all = frms_airpuff_all(:,any(frms_airpuff_all)); % remove the columns that contain zeros

    frms_airpuff_stay = [];
    frms_airpuff_run = [];
    for i = 1: length(airpuffon)
        if airpuffon(i)-before<1 || airpuffon(i)+after > length(speed)
            continue
        elseif sum(speed((airpuffon(i)-before):(airpuffon(i)-1))==0) == before
            frms_airpuff_stay = cat(1,frms_airpuff_stay, airpuffon(i)- before : airpuffon(i)+after);% trail*frames
        elseif sum(speed((airpuffon(i)-before):(airpuffon(i)-1))> 0) == before
            frms_airpuff_run = cat(1,frms_airpuff_run, airpuffon(i)- before : airpuffon(i)+after);% trial * frames
        end
    end
 % generate matrix for speeds 
    spd_airpuff_all = speed(frms_airpuff_all);
    spd_airpuff_all = spd_airpuff_all';
    
    spd_airTrig_run = speed(frms_airpuff_run);
    spd_airTrig_stay = speed(frms_airpuff_stay);
% averages and stes  
    spd_airTrigAve = mean(spd_airpuff_all);
    spd_airTrigAveSte = std(spd_airpuff_all)/sqrt(length(airpuffon));
    
    spd_airTrigAve_run = mean(spd_airTrig_run);
    spd_airTrigAveSte_run = std(spd_airTrig_run)/sqrt(size(spd_airTrig_run,1));
    
    spd_airTrigAve_stay = mean(spd_airTrig_stay);
    spd_airTrigAveSte_stay = std(spd_airTrig_stay)/sqrt(size(spd_airTrig_stay,1));
    
% figure    
    spd_fig = figure;
    x1 = 1: length(spd_airTrigAve); x2 = 1:length(spd_airTrigAve_run);
    subplot(1,3,1);
    shadedErrorBar(x1,spd_airTrigAve,spd_airTrigAveSte);
    ylabel ('speed');vline(period+1,'r');
    title(['airpuff trig average all trials',sessions{ii}]);
    subplot(1,3,2);
    shadedErrorBar(x2,spd_airTrigAve_run,spd_airTrigAveSte_run);
    xlabel('frames');vline(before+1,'r');
    title('running');
    subplot(1,3,3);
    shadedErrorBar(x2,spd_airTrigAve_stay,spd_airTrigAveSte_stay);
    vline(before+1,'r');
    title('stationary');
   
    saveas(spd_fig, [behav_dest '\' sessionID{ii} '_airpuffTrigAve']);
    save([behav_dest '\' sessionID{ii} '_behavAnalysis.mat' ],...
        'frms_airpuff_all','spd_airpuff_all','spd_airTrigAve','spd_airTrigAveSte',...
        'frms_airpuff_run','spd_airTrig_run','spd_airTrigAve_run','spd_airTrigAveSte_run',...
        'frms_airpuff_stay','spd_airTrig_stay','spd_airTrigAve_stay','spd_airTrigAveSte_stay','-append');
end
