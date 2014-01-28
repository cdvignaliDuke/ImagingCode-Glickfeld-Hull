function out = imAdaptHistEqJava(in,winsize)
%IMADAPTHISTEQJAVA
% OUT = IMADAPTHISTEQJAVA(IN,WINSIZE)
%
%   WRAPPER AROUND FILTER_RANKAK (MODIFICATION OF FILTER_RANK IJ PLUGIN)
%   DOES ADAPTIVE HISTOGRAM EQUALIZATION
%
%   CODE FROM AARON'S ADAPTHISTEQAK.M
%
%$ID

in = im2uint8(in);
si = size(in);
FR = Filter_RankAK;
FR.doIrescale = 1;
FR.doIrandomise = 0;
FR.doInewimage = 0;
FR.windowsize = winsize;

INIJ = ijarray2plus(in,'uint8');
FR.setup('',INIJ);
IP = INIJ.getProcessor;
disp('Processing adaptive histogram equalization in Java...please wait');
FR.run(IP);
disp('Adaptive histogram equalization complete');
pix = IP.getPixels;
out = subFixSign((reshape(pix,si(2),si(1)))');

%%%%%%%%%%%%%%%%

function out = subFixSign(in)
% maps -128 to 0 range of a signed integer to 128-255 of unsigned
in = double(in);
I = find(in<0);
in(I) = in(I)+256;
out = uint8(in);
