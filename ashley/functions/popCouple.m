function [pCpl,pCpl_FR] = popCouple(dff,frRateHz)
% population coupling = (1/stdev(f(cell i)) * integral[ f(cell i) * sum( f(all other cells) - (their mean f) ) ]
nc = size(dff,2);

pCpl = nan(1,nc);
pCpl_FRprod = nan(1,nc);
for i = 1:nc
    j = logical(ones(1,nc));
    j(i) = false;

    fi = dff(:,i);
    fj = dff(:,j);
    
    uj = mean(fj,1);
    fjresid = bsxfun(@minus,fj,uj);

    fjsum = sum(fjresid,2);
    
    fifj = fi.*fjsum;
    
    fifj_int = sum(fifj)/frRateHz;
    
    normfactor = 1/std(fi);
    
    ci = fifj_int/normfactor;

    pCpl(i) = ci;
    fi_fr = mean(fi);
    fj_fr = mean(sum(fj,2));
    pCpl_FRprod(i) = fi_fr.*fj_fr;
end
pCpl_FR = sqrt(abs(pCpl_FRprod));

end
