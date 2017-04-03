function [valTr, invTr] = matchTrialLengthIndex(val_trL,inv_trL,tars,minTrLMs,trType)
%%
if strcmp(trType,'short')
    val_ind = cellfun(@(x) find(x < minTrLMs), val_trL, 'unif',0);
    inv_ind = cellfun(@(x) find(x < minTrLMs), inv_trL, 'unif',0);
elseif strcmp(trType,'long')
    val_ind = cellfun(@(x) find(x > minTrLMs), val_trL, 'unif',0);
    inv_ind = cellfun(@(x) find(x > minTrLMs), inv_trL, 'unif',0);
elseif strcmp(trType,'all')
    val_ind = cellfun(@(x) find(x > 0), val_trL, 'unif',0);
    inv_ind = cellfun(@(x) find(x > 0), inv_trL, 'unif',0);    
end

if length(tars) == 2
    valTr = cell(1,length(tars));
    invTr = cell(1,length(tars));
    if any(cellfun(@(x) isempty(x), inv_ind)) 
        if ~all(cellfun(@(x) isempty(x), inv_ind)) 
            tar_ind = find(cellfun(@(x) isempty(x), inv_ind));
            if length(tar_ind) > 1
               valTr = cell(1,length(tars)); 
               invTr = cell(1,length(tars)); 
            else
                valTr{tar_ind} = [];
                invTr{tar_ind} = [];
                valTr{setdiff(1:2,tar_ind)} = val_ind{setdiff(1:2,tar_ind)};        
                invTr{setdiff(1:2,tar_ind)} = inv_ind{setdiff(1:2,tar_ind)};
            end
        end
    elseif any(cellfun(@(x) isempty(x), val_ind))
        tar_ind = find(cellfun(@(x) isempty(x), val_ind));
        if length(tar_ind) > 1
           valTr = cell(1,length(tars)); 
           invTr = cell(1,length(tars)); 
        else
            valTr{tar_ind} = [];
            invTr{tar_ind} = [];
            valTr{setdiff(1:2,tar_ind)} = val_ind{setdiff(1:2,tar_ind)};        
            invTr{setdiff(1:2,tar_ind)} = inv_ind{setdiff(1:2,tar_ind)};
        end
    else
        % match across target types if more than one
        nVal = cellfun(@(x) size(x,2), val_ind);
        nInv = cellfun(@(x) size(x,2), inv_ind);

        [nMatchVal, nMatchInv] = matchTrialRatio(nVal,nInv);

        for itar = 1:length(tars)
            nval_temp = nMatchVal(itar);
            ninv_temp = nMatchInv(itar);
            v_temp = val_ind{itar};
            inv_temp = inv_ind{itar};
            v_samp = randsample(nVal(itar),nval_temp);
            inv_samp = randsample(nInv(itar),ninv_temp);
            valTr{itar} = sort(v_temp(v_samp));
            invTr{itar} = sort(inv_temp(inv_samp));
        end
    end
else
    valTr{1} = val_ind;
    invTr{1} = inv_ind;
end
end