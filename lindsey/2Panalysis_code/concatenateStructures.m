function outS = concatenateStructures(data)
%% combining fields across structures

    fNames = fieldnames(data(1,1));
    nF = length(fNames);
    outS = struct([]);
    for iF = 1:nF
        fN = cell2mat(fNames(iF,:));
        for ic = 1:size(data,2)
            outS(ic).(fN) = data(1,ic).(fN);
            for iblock = 2:size(data,1)
                if size(outS(ic).(fN),1) == 1
                    outS(ic).(fN) = [outS(ic).(fN) data(iblock,ic).(fN)];
                else
                    outS(ic).(fN) = [outS(ic).(fN); data(iblock,ic).(fN)];
                end
            end
        end
    end