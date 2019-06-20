%% load data

clear;
sessions = '190430_img1023'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1023-190430_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave = load([image_analysis_dest 'getTC\' filename.name]);
tc_avg = TCave.tc_avg;

%load('190430_img1023_000_nPCA200_mu0.3_nIC80_thresh96_TCave.mat')

%%
% inputs:
% data
% threshold (number of std dev of the derivative)
% refractory period (in frames)
% min spike width (in frames)

% if want to use convolution
% tc_avg = conv2(tc_avg', ones(1,10)/10, 'same')';
% data, threshold, refractoryperiod, minspkwidth
spklogic = detectspikes(tc_avg, 2, 0, 0);
%the numbers right now give me the same thing as my code. if change the
%refractory period/min spike width as non zero, will get a much lower spike
%rate. 

%% plot spklogic
figure;
imagesc(spklogic)
xlabel('Cell #');
ylabel('Frame #');

%% plot example trace
c = 5;
figure;
subplot(3,1,1)
plot(tc_avg(:,c)); hold on;

spk2plot = tc_avg(:,c).*spklogic(:,c);
spk2plot(spk2plot==0) = NaN;

plot(spk2plot, 'ro', 'MarkerSize', 6)
xlim([1 300])

subplot(3,1,2)
plot(conv(tc_avg(:,c), ones(1,10)/10, 'same'))
xlim([1 300])

subplot(3,1,3)
plot(diff(tc_avg(:,c),1));
xlim([1 300])


%% calculate spike rates
% sampling period * number of frames
% subtract one frame for t=0

totaltime = (1/30)*(size(tc_avg,1)-1);

spkrates = sum(spklogic)/totaltime;

figure; 
scatter(1:length(spkrates), spkrates, 'bo')
ylabel('Spike Rate (Hz)');
xlabel('Cell #');


