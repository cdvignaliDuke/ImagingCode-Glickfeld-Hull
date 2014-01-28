function filtered = lowpass(timeCourses,Wn)
%LOWPASS
% FILTERED = LOWPASS(TIMECOURSES,WN)
% Wn is cutoff normalized frequency (default .3)

if nargin < 2
    Wn = .3;
end;

[b,a]=butter(5,Wn);

filtered = filtfilt(b,a,timeCourses);

return;

