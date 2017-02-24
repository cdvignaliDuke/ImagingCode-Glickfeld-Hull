function val_inv_cdf(valid,invalid,plotName)

h = cdfplot(valid);
h.Color = 'k';
hold on
h = cdfplot(invalid);
h.Color = 'c';
title(plotName);
ylabel('fraction of cells');
end