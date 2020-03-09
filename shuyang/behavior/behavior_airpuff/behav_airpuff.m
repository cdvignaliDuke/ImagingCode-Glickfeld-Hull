%This script is to analyze the RAW behavior data of running behavior during imaging
%The analysis will calculate average running speed, average duration of
%running periods, and average duration of stationary periods.
%This analysis will also graph speed across time, and make histograms of durations of each behavioral state.
%will save speed, cells during different behavioral states, median& std& ave of duration of each behavioral state.

%% Section I: set paths and create analysis folders for each session
%define the directory and files
clear;
folder = 'Z:\Data\behavior\airpuff\';
sessionID = '1042-191115';% this  variable name is confusing, this session ID is just tha date and the subject#, 
%there might be more than 1 sessions on a single subject on the same day
filename = dir([folder 'data-i' '*' sessionID  '*' ]);
for i = 1: size(filename,1)
    behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' sessionID '_' num2str(i)];
    if ~exist(behav_dest)
        mkdir(behav_dest);
    end
end

%% Section II : calculate for each frame,find relative frames for each behaviroal state and plot speed.
for i = 1:size(filename,1)
    file = [folder, filename(i).name];
    load(file);
    % later analysis will need to know if this session is a reverse session, and if yes, the frames of the reverse stimuli. 
    airpuffon = input.cTactileStimTurnedOn;
    airpuffon1 = airpuffon(~cellfun('isempty',airpuffon)); %it seems the airpuff sometimes is skipped, and the cell is empty so just extract the nonempty cells
    airpuffon1 = double(cell2mat(airpuffon1));
    lenairpuff = double(input.tactileStimulusDurationUs);
    lenairpuff = lenairpuff/1000000;
    lenframe = 33.3333;
    speed = calculate_speed_2P(input,lenframe);
    speed = double(speed);
    % delete the first frame because the first frame of imaging data is half black
    speed = speed(2:end);
    % find relative behavioral states and save to behavior analysis file
    smallestspd = ceil(1/lenframe*1000);% possible smallest speed: if the mice only moves 1unit during a frame, average speed of that frame is  (1/framelength) *1000
    [frames,frames_stay_cell, frames_bf_cell, frames_run_cell, frames_move_cell] = findFrames_behavStates_2P(speed,smallestspd);
    save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
        'airpuffon','speed','frames','frames_stay_cell','frames_bf_cell',...
        'frames_run_cell','frames_move_cell','airpuffon1');
    % plot speed and save figure
    fig_speedtc = figure;
    plot(frames, speed); hold on;
    if isempty(airpuffon1) == 0
    vline(airpuffon1, 'r');
    end
    title (['average speed every frame(33ms)', '  ', sessionID '-' num2str(i)]);
    xlabel ('frames');
    ylabel ('speed(pulses/s)');
    text(10,-20,['airpuff length ' num2str(lenairpuff) 's']);
    saveas(fig_speedtc,[behav_dest '\' sessionID '_' num2str(i) '_speed']);
end
% check if found frames are what I want to find
figure;plot(speed);
hold on;
plot(cell2mat(frames_bf_cell), smallestspd*ones(1,length(cell2mat(frames_bf_cell))),'r.');


%% SECTION III: draw distribution of running duration
%num_frames_run_all = {};
num_frames_run = [];
imgfreq = 30; %!!!!!!!!!!!!!!!!!!!!!!make sure it's the right frequency!
for i =  1:size(filename,1)
    for n = 1:size(frames_run_cell,2)
    num_frames_run = [num_frames_run size(frames_run_cell{n},2)];
    end
    aveDuration_runFrames = mean(num_frames_run);
    medianDura_runFrames = median(num_frames_run);
    stdDura_runFrames = std(num_frames_run);
    % frames to seconds
    aveDura_run = aveDuration_runFrames/imgfreq;
    medianDura_run = medianDura_runFrames/imgfreq;
    stdDura_run = stdDura_runFrames/imgfreq;
    
    save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
        'aveDuration_runFrames', 'medianDura_runFrames','stdDura_runFrames', ...
        'aveDura_run','medianDura_run','stdDura_run','-append');
    % plot distribution of all running windows' duration
    fig_runDuradist = figure;
    edges = (0:1:max(num_frames_run/imgfreq));
    histogram(num_frames_run/imgfreq,edges);
    title (['distribution of running duration', '  ',filename(i).name]);
    xlabel ('time length (seconds)');
    ylabel ('frequency');
    vline(medianDura_runFrames/imgfreq, 'r');
    saveas(fig_runDuradist,[behav_dest '\' sessionID '_' num2str(i) '_runDuraDist']);
end


%% Section IV: Calculate average stationary duration:
num_frames_stay = [];
imgfreq = 30; %!!!!!!!!!!!!!!!!!!!!!!make sure it's the right frequency!

for i = 1:size(filename,1)
    for k = 1: size(frames_stay_cell,2)
        num_frames_stay = [num_frames_stay size(frames_stay_cell{k},2)];
    end
    aveDuration_stayFrames = mean(num_frames_stay);
    medianDura_stayFrames = median(num_frames_stay);
    stdDura_stayFrames = std(num_frames_stay);
    % frames to seconds
    aveDura_stay = aveDuration_stayFrames/imgfreq;
    medianDura_stay = medianDura_stayFrames/imgfreq;
    stdDura_stay = stdDura_stayFrames/imgfreq;
    
    save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
        'aveDuration_stayFrames', 'medianDura_stayFrames','stdDura_stayFrames',...
        'aveDura_stay','medianDura_stay','stdDura_stay','-append');
    %plot distribution of number of frames in each stationary state
    fig_staydist = figure;
    edges = (0:1:max(num_frames_stay/imgfreq));
    histogram(num_frames_stay/imgfreq,edges);
    title (['distribution of stationary duration', '  ',filename(i).name]);
    xlabel ('time length (seconds)');
    xlim([1 500]);
    ylabel ('frequency');
    vline(medianDura_stay, 'r');
    saveas(fig_staydist,[behav_dest '\' sessionID '_' num2str(i) '_stayDistribution' ]);
end


%% Section V: calculate average speed and ave&median duration of running above the runThreshold(mean)
% calculate the average speed to get running threshold

imgfreq = 30; %!!!!!!!!!!!!!!!!!!!!!!make sure it's the right frequency!
for i = 1:size(filename,1)
    frames_run = frames(cell2mat(frames_run_cell));
    speed_run = speed(frames_run);
    aveRunSpeed = mean(speed_run);
    %threshold of reverse stimuli.
    runThreshold = roundn(aveRunSpeed,1); %round averun to the nearest multiple of 10.Need to do this because all speed values are multiple of 10, so need a multiple of 10 for setting the
   
    % find above threshold running durations
    dura_runAbvThre = [];
    for p = 1: size(frames_run_cell,2)
        dura_runAbvThre(p) = sum((speed(frames_run_cell{p}) >= runThreshold));
    end
    % frames to seconds
    aveDura_runAbvThre = mean(dura_runAbvThre)/imgfreq;
    medianDura_runAbvThre = median(dura_runAbvThre)/imgfreq;
    stdDura_runAbvThre = std(dura_runAbvThre)/imgfreq;
    save([behav_dest '\' sessionID '_' num2str(i) '_behavAnalysis.mat' ],...
        'runThreshold', 'aveDura_runAbvThre','medianDura_runAbvThre','stdDura_runAbvThre', '-append');
    %plot distribution of number of frames in each stationary state
    fig_ABTspeeddist = figure;
    edges = (0:1:max(dura_runAbvThre/imgfreq));
    histogram(dura_runAbvThre/imgfreq,edges);
    title (['distribution of abvThreRun duration', '  ',filename(i).name]);
    xlabel ('time length (seconds)');
    ylabel ('frequency');
    vline(medianDura_runAbvThre, 'r');
    saveas(fig_ABTspeeddist,[behav_dest '\' sessionID '_' num2str(i) '_ATR Distribution' ]);
end

