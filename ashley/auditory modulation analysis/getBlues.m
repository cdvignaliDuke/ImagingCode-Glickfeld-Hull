function continuousBlueColorSet = getBlues(nType)
    colorset = brewermap(nType-1,'*Blues');
    delayColors = colorset(1:end-1,:);
    continuousBlueColorSet = cat(1,delayColors,[0.5 0.5 0.5],[0 0 0]);
end