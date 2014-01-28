function figH = stack_newplot
%
%$Id: stack_newplot.m 97 2008-03-17 20:14:13Z histed $

figH = findobj(0, 'Tag', 'MHStackPlotFig');
if isempty(figH)
    figH = figure('Name', 'Stack plotter', ...
                  'Tag', 'MHStackPlotFig');
else
    clf;
end

