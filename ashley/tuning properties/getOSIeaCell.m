function osi = getOSIeaCell(oriFit,orientations)
fitOrthos = circshift(orientations,length(orientations)/2);

[prefResp, prefInd] = max(oriFit,[],1);
nc = size(oriFit,2);
osi = nan(1,nc);
for i = 1:nc
    pref = prefResp(i);
    orthInd = fitOrthos(prefInd(i));
    orth = oriFit(orthInd,i);
    osi(i) = (pref-orth)./(pref+orth);
end

end