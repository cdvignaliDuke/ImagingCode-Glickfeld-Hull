function [inc, dec, none, grp] = auROCsignifGroups(auroc,rst)

signif = find(rst);
inc = intersect(find(auroc > 0.5),signif); 
dec = intersect(find(auroc < 0.5),signif); 
none = setdiff(1:length(auroc),[inc,dec]);

grp = {inc;dec;none};

end
