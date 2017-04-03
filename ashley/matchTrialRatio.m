function [nMatchVal, nMatchInv] = matchTrialRatio(nVal,nInv)

invRatio = nInv(1)/nInv(2);
valRatio = nVal(1)/nVal(2);

nMatchVal = zeros(1,2);
nMatchInv = zeros(1,2);
gt_ind = gt(nInv,nVal);
if all(gt_ind)% more invalid trials for each target type    
    if invRatio > valRatio
        nMatchInv(1) = nInv(1);
        nMatchInv(2) = round(nInv(1) * (1/valRatio));
    else
       nMatchInv(2) = nInv(2);
       nMatchInv(1) = floor(nInv(2) * valRatio);
    end  
    nMatchVal = nVal;    
elseif gt_ind(2) % more invalid trials for last target type
    if invRatio > valRatio
        nMatchVal(1) = nVal(1);
        nMatchVal(2) = round(nVal(1) * (1/invRatio));
    else
       nMatchVal(2) = nVal(2);
       nMatchVal(1) = floor(nVal(2) * invRatio);
    end  
    nMatchInv = nInv;
elseif ~any(gt_ind) || gt_ind(1)% more valid trials for each target type OR 
%                             more invalid trials for first target type
    if invRatio > valRatio
        nMatchVal(1) = nVal(1);
        nMatchVal(2) = round(nVal(1) * (1/invRatio));
    else
       nMatchVal(2) = nVal(2);
       nMatchVal(1) = floor(nVal(2) * invRatio);
    end
    nMatchInv = nInv;
end
end