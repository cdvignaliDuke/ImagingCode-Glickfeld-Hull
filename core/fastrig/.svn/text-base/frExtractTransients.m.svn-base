function t = frExtractTransients(timeCourses,varargin)
% T = FREXTRACTTRANSIENTS(TIMECOURSES,OPTIONS)

%% default arguments
defaultopts = {'fs',frGetFrameRate,'tau',0.5,'polyN',4,'frameLen', 21,'threshold',5};

options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);

for iarg = 1:2:length(varargin)
    options.(varargin{iarg}) = varargin{iarg+1};
end

t.timeCourses = timeCourses;
[t.nsamples,t.ncells]=size(timeCourses);
t.options = options;

if nargin == 0
    disp(options);
end
    
%% match filtering, removes some noise
lambda = 1/options.tau;
t.ts = [0:64]/options.fs; % in seconds at options.fs hz
t.filt = exppdf(t.ts,options.tau)/options.fs; % normalizes sum 1 to
t.matched = filter(t.filt,1,timeCourses);

%% Savitzky-Golay polynomial filter, remove noise and differentiate
[t.interp,t.diff]=sgolaydiff(t.matched,options.polyN,options.frameLen);

%% extract transients
outliers = find(abs(t.diff)>0.5*repmat(iqr(t.diff),t.nsamples,1) | t.diff > 0.0);
temp = t.diff;temp(outliers)=nan;
lessbiased = 2*nanstd(temp);

t.thresholds = options.threshold * repmat(lessbiased,t.nsamples,1);
t.bursts = t.diff>t.thresholds; % true during bursts
t.onsets = diff(t.bursts)==1;
t.offsets = diff(t.bursts)==-1;

onsetlist = find(t.onsets);
offsetlist = find(t.offsets);

t.bursttimes = zeros(size(t.bursts));
t.amplitudes = zeros(size(t.bursts));

for index = 1:length(onsetlist)
%    [onsetlist(index),offsetlist(index)]  
    if ~isempty(onsetlist(index):offsetlist(index));
        t.amplitudes(onsetlist(index))=max(t.diff(onsetlist(index):offsetlist(index)));
        t.bursttimes(onsetlist(index))=1;
    end
end

return;

%%
figure;clf;
tcOffsetPlot(1:t.nsamples,t.diff(:,1:8:end),3);
hold on
tcOffsetPlot(1:t.nsamples,t.thresholds(:,1:8:end),3);
xlim([0 5000])


%%
figure;
tcOffsetPlot(1:t.nsamples,t.diff(:,1:8:end),3);
hold on
tcOffsetPlot(1:t.nsamples,t.diff(:,1:8:end),3);
xlim([0 5000])
