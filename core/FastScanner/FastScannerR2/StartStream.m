function p=StartStream(object, event)
% Starts streaming timer. Called upon acquisition trigger

global fs;

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

fs.DAQ.Streaming = 1;

tic;
start(fs.StreamTimer);

