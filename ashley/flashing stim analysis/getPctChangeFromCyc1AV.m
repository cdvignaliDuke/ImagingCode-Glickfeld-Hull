function pctChangeFromCyc1AV = getPctChangeFromCyc1AV(avgRespEaCycAV,...
    normCondition)

normValue = avgRespEaCycAV(1,normCondition);

pctChangeFromCyc1AV = avgRespEaCycAV./normValue;

end