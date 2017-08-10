function respAxLim = minMaxDFFRespAxLim(allRespVals,buffer)
    resp1 = min(allRespVals)-buffer;
    resp2 = max(allRespVals)+buffer;
    respAxLim = [resp1 resp2];
end