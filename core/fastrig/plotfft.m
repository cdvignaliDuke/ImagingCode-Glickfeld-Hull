function h = plotfft(xx,nfft,Fs,varargin);
%PLOTFFT
% [YY,FF] = PLOTFFT(XX,NFFT,FS);

if nargin < 2 | nargin > 1 & isempty(nfft)
    nfft = length(xx)
end

if nargin < 3
    Fs = frGetFrameRate;
end
   
if isempty(nfft)
    nfft = length(xx)
end;

nseg = floor(length(xx)/nfft);

xx = reshape(xx(1:nseg*nfft), [nfft nseg]);

% frequency vector 
ff = Fs/2*linspace(0,1,nfft/2);

XX = fft(xx,nfft)/sqrt(nfft);

% amplitude spectrum
yy = mean(sqrt(XX.*conj(XX)),2);

% power in dB
% db = 10*log(yy);

hh=plot(ff,yy(1:nfft/2),varargin{:});
xlabel('Frequency (Hz)');

yl = ylim;

if nargout > 0
    h = hh;
end;

return;
