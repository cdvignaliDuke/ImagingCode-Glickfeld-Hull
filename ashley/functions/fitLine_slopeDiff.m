function [pVal,slopes] = fitLine_slopeDiff(x1,y1,x2,y2);
    
n1 = length(x1);
n2 = length(x2);
grpID = cat(1,ones(n1,1),ones(n2,1).*2);
x = cat(1,x1,x2);
y = cat(1,y1,y2);
tbl = table(x,y,grpID);
fit = fitlm(tbl,'x~y*grpID');
tbl_anova = anova(fit);
pVal = tbl_anova.pValue('y:grpID');

end