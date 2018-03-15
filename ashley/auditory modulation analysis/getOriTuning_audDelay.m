function [vonMisesFit,fitReliability,R_square] = ...
    getOriTuning(responses,responses_resamp,orientations,nBoot)

[vonMisesFit,~,fitReliability,R_square] = vonmisesReliableFit(responses,...
    responses_resamp,orientations,nBoot);

end