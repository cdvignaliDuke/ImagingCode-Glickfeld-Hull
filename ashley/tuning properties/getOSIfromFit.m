function OSI = getOSIfromFit(oriFit)
    
    
    [maxResp,oriPref] = max(oriFit,[],1);
    nullOris = circshift(0:180,90);
    nc = length(oriPref);
    nullResp = nan(1,nc);
    for icell = 1:nc
        nullResp(icell) = oriFit(nullOris == oriPref(icell),icell);
    end
    OSI = (maxResp-nullResp)./(maxResp+nullResp);

end