function hh = tcOffsetPlot(tt,timeCourses,delta,varargin)
%TCOFFSETPLOT Plot time courses with offset
% h = tcOffsetPlot(timeCourses) plots time courses
% h = tcOffsetPlot(tt,timeCourses) specifies x axes. Default is sample #.
% h = tcOffsetPlot(tt,timeCourses,delta) specifies offset between traces.
%  Default is three times the average standard deviation of the traces.
% h = tcOffsetPlot(tt,timeCourses,delta, ...) additional arguments to plot

if nargin < 2
    timeCourses = tt;
    [nsamples,ncells]=size(timeCourses);
    tt = [0:nsamples-1];
end

if nargin < 3 | isempty(delta)
    if length(timeCourses)==1
        delta = timeCourses;
        timeCourses = tt;
        [nsamples,ncells]=size(timeCourses);
        tt = [0:nsamples-1];
    else
        delta = double(5*nanmean(nanstd(timeCourses)));
    end
end

[nsamples,ncells]=size(timeCourses);

% timeCoursesMinusDC = tcRemoveDC(timeCourses);

offsets = repmat(0:delta:(ncells-1)*delta,nsamples,1);

h =plot(tt,timeCourses+offsets,varargin{:});

if ncells >1 
    ylim([-delta delta*(ncells)]);
end;

set(h(2:end),'handlevisibility','off');

box off;
% if nargout > 1     hh = h;end
hh = h;
return;