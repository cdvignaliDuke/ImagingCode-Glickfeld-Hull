function [normedResponsesPrimary, normedResponsesSecondary] = ...
    norm2PeakOfOtherResp(responsesOfPrimary, responsesOfSecondary)
% responsesOfPrimary and responsesOfSecondary are matrices of the same
% size, observations x neurons

primaryPeak = max(responsesOfPrimary,[],1);
normedResponsesPrimary = bsxfun(@rdivide, responsesOfPrimary, primaryPeak);
normedResponsesSecondary = bsxfun(@rdivide, responsesOfSecondary, primaryPeak);

end