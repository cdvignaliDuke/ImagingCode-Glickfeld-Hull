function saveslides(fn,list)
%SAVESLIDES
% SAVESLIDES(FILENAME,LIST) where FILENAME is the name of the ppt
% presentation and LIST is a list of figure handles. Each figure pasted
% onto one slide.
%
% See saveppt2

% Wrapper needed because saveppt2 only allows to append one slide per call.

% if exist(fn)
%     error('file already exists');
% end

notes = sprintf('Printed by %s on %s. -- %s -- %s',...
                username,datestr(now),hostname,pwd)

for ind = 1:length(list)
%    name= get(list(ind),'name');
    fprintf('Saving Figure %i\n',list(ind));
    figure(list(ind));
    saveppt2(fn,'figure',list(ind),...
            'stretch','off','device','bitmap','scale','on','notes',notes);
end

return