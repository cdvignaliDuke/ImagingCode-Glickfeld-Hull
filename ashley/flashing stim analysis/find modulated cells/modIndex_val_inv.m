% modulation index = (val-inv)/(val+inv)
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
%hits only
val = resp_hvsfa_all;
inv = resp_favsh_all;

mi_hvsfa_all = (val-inv)./(val+inv);


%all
val = resp_val_all;
inv = resp_inv_all;

mi_valvsinv_all = (val-inv)./(val+inv);

cdfMI = figure; setFigParams4Print('landscape')
suptitle('attention modulation index')
subplot(1,2,1)
f = cdfplot(mi_hvsfa_all);
f.Color = color;
f.LineWidth = 3;
xlim([-20 20])
title(['H-v vs H-inv'])
xlabel('(val-inv)/(val+inv)')
ylabel('fraction of cells')
axis square

subplot(1,2,2)
f = cdfplot(mi_valvsinv_all);
f.Color = color;
f.LineWidth = 3;
xlim([-20 20])
title(['all-val vs all-inv'])
xlabel('(val-inv)/(val+inv)')
ylabel('fraction of cells')
axis square


print([fnout 'catch_align_miCDF' datasetStr '.pdf'], '-dpdf','-fillpage')
