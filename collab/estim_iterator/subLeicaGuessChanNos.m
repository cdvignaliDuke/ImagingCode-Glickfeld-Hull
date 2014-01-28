function [greenChanNo redChanNo] ...
    = subLeicaGuessChanNos(dataPath, seriesName)
%
%$Id: subLeicaGuessChanNos.m 394 2008-11-21 07:02:34Z histed $

list = dir(dataPath);
filelist = {list.name};
if isempty(filelist)
    error('No tif files found at %s', dataPath);
end

% pull out the list of channel numbers
tRE = sprintf('%s_t[0-9]+_ch([0-9]+)\\.tif+', ...
              seriesName);
tokC = regexp(filelist, tRE, 'tokens');
r = cat(1,tokC{:});
channelList = str2double(cat(1,r{:}));
unqChans = unique(channelList);

% check if there are one or two channels and select one
if length(unqChans) == 1
    assert(unqChans == 0, 'bug?');
    redChanNo = [];
    greenChanNo = 0;
else
    assert(all(unqChans(:) == [0 1]'));
    assert(length(channelList)/2==sum(channelList==0), 'bug?');
    redChanNo = 0;
    greenChanNo = 1;
end
