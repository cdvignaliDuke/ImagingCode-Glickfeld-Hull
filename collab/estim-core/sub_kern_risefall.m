function [outKern] = sub_kern_risefall(sampTimeMs, fallTauMs, riseTauMs)
%$Id: sub_kern_risefall.m 205 2008-05-02 03:37:09Z histed $

%%% make kernel

%sampTimeMs = lineTimeMs;
%fallTauMs = 125;
%riseDFOF = 0.04;
%riseTauMs = 4;
xPtsSamp = 0:(fallTauMs*8./sampTimeMs);
kernS = -exp(-xPtsSamp./(riseTauMs/sampTimeMs)) + exp(-xPtsSamp./(fallTauMs/sampTimeMs));
outKern = kernS ./ max(kernS);  % norm to max

