function AVcdf(visual,auditory,plotName)

h = cdfplot(visual);
h.Color = 'g';
hold on
h = cdfplot(auditory);
h.Color = 'k';
title(plotName);
ylabel('fraction of cells');
end