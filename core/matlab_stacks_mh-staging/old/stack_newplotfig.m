function figH = stack_newplotfig
%
%$Id: stack_newplotfig.m 97 2008-03-17 20:14:13Z histed $

global MHSTACK;

figOpts = {};
if isfield(MHSTACK, 'DefaultFigurePosition')
    figOpts = {figOpts{:}, 'Position', MHSTACK.DefaultFigurePosition};
end

figH = figure(figOpts{:});

% $$$ figH = findobj(0, 'Tag', 'MHStackPlotFig');
% $$$ if isempty(figH)
% $$$     figH = figure('Name', 'Stack plotter', ...
% $$$                   'Tag', 'MHStackPlotFig');
% $$$ else
% $$$     clf;
% $$$ end

