% modulation index = (vis-aud)/(vis+aud)
if cellsInd == 13 | cellsInd == 14
    color = [0 0 0];
elseif cellsInd == 2
   color = [0.5 0.5 0.5];
elseif cellsInd == 3
   color = 'r';
elseif cellsInd == 4
   color = 'g';
elseif cellsInd == 5
   color = 'b';
end
%early
vis = earlyResp_vis;
aud = earlyResp_aud;

mi_early = (vis-aud)./(vis+aud);


%late
vis = lateResp_vis;
aud = lateResp_aud;

mi_late = (vis-aud)./(vis+aud);

subplot(5,2,9)
f = cdfplot(mi_early);
f.Color = color;
f.LineWidth = 3;
xlim([-20 20])
title(['early mi'])
xlabel('(vis-aud)/(vis+aud)')
ylabel('fraction of cells')
axis square

subplot(5,2,10)
f = cdfplot(mi_late);
f.Color = color;
f.LineWidth = 3;
xlim([-20 20])
title(['late mi'])
xlabel('(vis-aud)/(vis+aud)')
ylabel('fraction of cells')
axis square


print([fnout 'catch_align_miCDF' datasetStr '.pdf'], '-dpdf','-fillpage')
