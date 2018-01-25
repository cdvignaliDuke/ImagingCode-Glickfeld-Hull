function DSI = getDSIfromFit(dirFit)
    
    
    [maxResp,dirPref] = max(dirFit,[],1);
    nullDirs = circshift(0:360,180);
    nc = length(dirPref);
    nullResp = nan(1,nc);
    for icell = 1:nc
        nullResp(icell) = dirFit(nullDirs == dirPref(icell),icell);
    end
    DSI = (maxResp-nullResp)./(maxResp+nullResp);

end