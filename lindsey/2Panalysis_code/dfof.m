function dfof_out = dfof(data,f);
    % data - raw tc
    % f - average
    dfof_out = bsxfun(@rdivide, bsxfun(@minus, data, f), f);
end