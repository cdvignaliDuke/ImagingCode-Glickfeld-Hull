function y = tcLowCut (x, cut_off_period, filter_type, boundary)
%TCLOWCUT Removes low frequencies while keeping DC
%   Y = TCLOWCUT (X, CUT_OFF_PERIOD, FILTER_TYPE, BOUNDARY)
%
% low cut filter for time courses
% DC will not be subtracted. 
% filter_type: 
%   if this arg is not specified, or specified as '',
%   filter cuts sharply at the cut_off period,
%   so ripples may happen.
%   if this is specified as 'gaussian',
%   this filter cuts with a gaussian window,
%   so ripples will be suppressed.
% boundary:
%   if this is not specified or set as 0, 
%   some artifacts can happen around the boundaries,
%   when the first and last values x the data are very different,
%   because FFT assumes periodic boundary condition.
%   if this is set as 1, this problem will be managed
%   by making a time-course as periodic. 
%
% examples
%   tcLowCut (x, cut_off_period)
%       low-cut sharply at cut_off_period
%   tcLowCut (x, cut_off_period, 'gaussian')
%       low-cut by gaussian around cut_off_period
%   tcLowCut (x, cut_off_period, '', 1)
%       low-cut sharply at cut_off_period, considering boundaries
%   tcLowCut (x, cut_off_period, 'gaussian', 1)
%       low-cut by gaussian around cut_off_period, considering boundaries
%
%
%   Kenichi Ohki  2004. 9. 8. 

if nargin <= 2
    filter_type = '';
end

if nargin <= 3
    boundary =0;
end

[m,n]=size(x);

if (n>1) && (m>1)
    y = x;
    for i=1:n  % loop over columns
       y(:,i) = tcLowCut(x(:,i),cut_off_period, filter_type, boundary);
    end
    return
end
    
N = length(x);
filter_size = N/cut_off_period;

if boundary == 1
    N=N*2;
    filter_size = filter_size *2;
    x=[x; flipud(x)];
end

switch lower(filter_type)
    case ('gaussian')
        filter=fftshift(1-gausswin(N+1, N/(filter_size*2)));
        filter=[1; filter(2:N)];
    case ('nodc');
        filter=ones(N,1);
        filter_size=round(filter_size);
        filter([2:2+filter_size-1,N-filter_size+1:N])=0;
        filter(1)=0;
    otherwise
        filter=ones(N,1);
        filter_size=round(filter_size);
        filter([2:2+filter_size-1,N-filter_size+1:N])=0;
end

k=fft(x);
k=k.*filter;

y=real(ifft(k));

if boundary == 1
    y = y(1:N/2);
end
