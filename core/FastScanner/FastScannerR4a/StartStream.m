function p=StartStream(object, event)
% Starts streaming timer. Called upon acquisition trigger

global fs;

fs.DAQ.aoCurrentTriggerValue = fs.DAQ.aoAcqTriggerStreamValue;

%fs.timing.startime = cputime;
fs.timing.startime = clock;

% delete any left-over timer
timers = timerfind;
for itimer = 1:length(timers);
    if strcmp(func2str(timers(itimer).TimerFcn),'ProcessStream');
        stop(timers(itimer));
        delete(timers(itimer));
    end
end

fs.StreamTimer = timer('TimerFcn',@ProcessStream, 'Period', 5e-3,'ExecutionMode','FixedRate');

fs.sAcquired=zeros(100000,1);
fs.sAvailable=zeros(100000,1);
fs.totaltime=zeros(100000,1);
fs.tictoc=zeros(100000,1);
fs.sAvailableFrames=zeros(100000,1);

fs.DAQ.firstdone=0;
fs.ph=[];
fs.period=[];

fs.map = fopen(fs.DAQ.MemMapFile,'W');

fs.DAQ.Streaming = 1;
tic;
start(fs.StreamTimer);

return;
