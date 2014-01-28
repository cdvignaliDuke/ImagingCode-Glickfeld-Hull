function drawSpectrogram(varargin)
%DRAWSPECTROGRAM
% Computes and draw spectrogram as image. This is needed because SPECTROGRAM
% uses SURF to draw
% SPECTROGRAM(X)
% S = SPECTROGRAM(X,WINDOW)
% S = SPECTROGRAM(X,WINDOW,NOVERLAP)
% S = SPECTROGRAM(X,WINDOW,NOVERLAP,NFFT)
% S = SPECTROGRAM(X,WINDOW,NOVERLAP,NFFT,Fs)
%
% see SPECTROGRAM for help.

[s,f,t,p]=spectrogram(varargin{:});
imagesc(f,t,10*log10(p'));
box off;

set(gca,'ydir','normal')
xlabel('Frequency (Hz)');
ylabel('Time (s)');

return;
