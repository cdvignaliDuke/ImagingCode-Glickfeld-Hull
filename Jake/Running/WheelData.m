% Re adapted for 9/9 data   
% running wheel script
% shows distance run for all time periods
% this code assumes the duration of ITI is equal to the duration of one of
% the spCounters
% HARDCODED for an ISI = (1/2) the duration of each spCounter. 

%% Section1 extract running wheel data to a vector and graph
mworks = 'Z:\Data\WidefieldImaging\GCaMP\140909_img14_1\data-img14-140909-1107';
fname = 'img14_140909';  %for output only
pname = 'Z:\Data\WidefieldImaging\GCaMP\140909_img14_1';
load(mworks);
ntrial = size(input.spCounter1, 2);   %WARNING: sometimes it will collect an extra trial. If this happens your analysis is fucked. Try setting ntrial = input.stopAfterNTrials

RunData = NaN(1,[ntrial*(input.nScansOn+input.nScansOff+(input.postSoundPauseDurationMs/input.frameImagingRateMs))]); %creates empty vector for data. Final destination of running wheel data
IntData = NaN(ntrial ,input.nScansOn+input.nScansOff+(input.postSoundPauseDurationMs/input.frameImagingRateMs)); %creates empty matrix to hold wheel pulse data

IntData(:,1) = cell2mat(input.tITIWheelCounter)'; %WARNING: These four lines are hardcoded for 9/9 data. Rewrite for new versions of cbstim
IntData(:,2) = cell2mat(input.tITIWheelCounter)'; %adds ITI and ISI pulses to the beginning of the intdata matrix
IntData(:,3) = cell2mat(input.tISIWheelCounter)';
IntData(:,4) = cell2mat(input.tISIWheelCounter)';

aa=1;
for b = 5:2:size(IntData,2);     %for loop to load wheel data into intdata from spCounters WARNING: starting valuee of b is hard coded for 9/9 data. rewrite. 
    IntData(:,b) = cell2mat_padded(eval(['input.spCounter' num2str(aa)]));
    IntData(:,b+1) = cell2mat_padded(eval(['input.spCounter' num2str(aa)]));
    aa = aa +1;
    %strcat would be better to do this sort of thing... or sprintf and %s
end

for d = 1:ntrial;
    RunData([1+(d-1)*size(IntData,2)]:[size(IntData,2)*d]) = IntData(d,:);
end

%generates a vector to serve as the X axis (time) for plotting the data. Each unit will equal the framerate (in Ms) times the # of frames in each of the 10 spCounters.
indmat = NaN(1,size(RunData,2));   
for a =  1:[ntrial*size(IntData,2)];
    indmat(a) = a;
end
plot(indmat, RunData)
clear a;

%% Section2 cut away frames from nScansOff sections 
%HARDCODED LINES IN HERE 500
imdata = readtiff(pname);
imdata2 = NaN(size(imdata, 1), size(imdata, 2), [ntrial*(input.nScansOn+(input.postSoundPauseDurationMs/500))]);  
for a = 1:ntrial;
%SHOULD BE SUBTRACTING NSCANON NOT NSCANOFF
    imdata2(:,:,[1+(a-1)*(input.nScansOn+input.postSoundPauseDurationMs/500):(a*(input.nScansOn+input.postSoundPauseDurationMs/500))])=imdata(:,:,[(a*(input.nScansOn+input.nScansOff+input.postSoundPauseDurationMs/500))-input.nScansOff:(a*(input.nScansOn+input.nScansOff+input.postSoundPauseDurationMs/500))]);
end

size(imdata2)
[row, col] = find(isnan(imdata2));
row
col     %both should be empty. If not then the rows and columns designated will reveal the location of any pesky NaNs

%% Section3 Set threshold value for running/not-running states and use runData to split imdata2 into seperate matrices by running state 
%for this experiment (9/4/14) I have defined running at >50 pulses per
%500ms frame. Not-running as <30 pulses per 500ms frame. This animal seemed
%to run a lot and never rest so this may be more of a running-fast vs
%running-slow division. 
clear a b c;
b=1;
c=1;
aa=0;
bb=0;
for a = 1:size(RunData,2);
    if RunData(:,a) > 50;
        aa=aa+1;
    elseif RunData(:,a) < 30;
        bb=bb+1;
    end
end
runIms = NaN(size(imdata2,1),size(imdata2,2),aa);
nonRunIms = NaN(size(imdata2,1),size(imdata2,2),bb);
for a = 1:size(RunData,2);
    if RunData(:,a) > 50; 
        runIms(:,:,b)=imdata2(:,:,a);
        b=b+1;
    elseif RunData(:,a) < 30;
        nonRunIms(:,:,c)=imdata2(:,:,a);
        c=c+1;
    end
end
runIms = mean(runIms,3);
nonRunIms = mean(nonRunIms,3);

%% Section4 write tifss
writetiff(runIms, [pname '\RunImsTiff']); 
writetiff(nonRunIms, [pname '\nonRunImsTiff']);







