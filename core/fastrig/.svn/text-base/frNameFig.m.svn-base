function str = frNameFig(fig,expt,label)

if nargin < 3
    label = expt;
    expt = fig;
    fig = gcf;
end

str = sprintf('%s-%s-%s',expt.animal,expt.name,label);

set(gcf,'numbertitle','off')
set(gcf,'name',str);

return;
