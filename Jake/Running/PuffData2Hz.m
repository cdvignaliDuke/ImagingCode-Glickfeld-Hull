% Meant to collect data from puff trials
% HARDCODED for an ISI = (1/2) the duration of each spCounter. 

%% Section1 
mworks = 'Z:\Data\WidefieldImaging\GCaMP\ToneAndPuff\140918_img10_1\data-img10-140918-1424';
fname = 'img10_140918';  %for output only
pname = 'Z:\Data\WidefieldImaging\GCaMP\ToneAndPuff\140918_img10_1';
subjname = 'img10';
date = '140918';
load(mworks);

ntrial = size(input.spCounter1, 2) -1;   %WARNING: sometimes it will collect an extra trial. If this happens your analysis is fucked. Try setting ntrial = input.stopAfterNTrials
nFramesTrial = (input.nScansOn + input.nScansOff + input.preSoundPauseDurationMs + input.postSoundPauseDurationMs)/input.frameImagingRateMs;

%% Section2 cut away frames from nScansOff sections 
%HARDCODED LINES IN HERE 500
imdata = readtiff(pname);
imBaseline = NaN(size(imdata,1), size(imdata,2), ntrial);
imTone = NaN(size(imdata,1), size(imdata,2),ntrial);
imPuff = NaN(size(imdata,1), size(imdata,2),ntrial);

for a = 1:ntrial;
    imBaseline(:,:,a)=imdata(:,:,[1+((a-1)*ntrial)]);
    imTone(:,:,a)=imdata(:,:,[22 + (a-1)*ntrial]);      
    imPuff(:,:,a)=imdata(:,:,[23 + (a-1)*ntrial]);
end

imBaseline=mean(imBaseline,3);
imTone=mean(imTone,3);
imPuff=mean(imPuff,3);



%% Section4 write tifss
writetiff(imBaseline, [pname '\' subjname 'BaselineTiff' date]); 
writetiff(imTone, [pname '\' subjname 'ToneTiff' date]);
writetiff(imPuff, [pname '\' subjname 'PuffTiff' date]);






