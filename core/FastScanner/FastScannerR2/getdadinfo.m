function [ndads,nframes]=getdadinfo(basepath)

ndads = length(dir([basepath '*']));
nframes = (ndads-1)*150;

return;

