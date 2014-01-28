function [modelout] = addpiecewiseeventsSB(model)
% addpiecewiseeventsSB: In order to switch pieceswiseSB constructs at
% precisely the right instants the trigger functions of all piecewise
% statements in a model are additionally implemented as events with dummy
% assignments.
%
% USAGE:
% ======
% [modelout] = addpiecewiseeventsSB(model)       
%
% model: SBmodel

% Information:
% ============
% Copyright (C) 2008 Henning Schmidt, henning@sbtoolbox2.org
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
% USA.

if ~isSBmodel(model),
    error('Function only defined for SBmodels.');
end
ms = struct(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND ALL PIECEWISE TRIGGERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triggers = {};
% Check in ODEs
for k=1:length(ms.states),
    addtriggers = gettrigger(ms.states(k).ODE);
    triggers = {triggers{:}, addtriggers{:}};
end
% Check in variables
for k=1:length(ms.variables),
    addtriggers = gettrigger(ms.variables(k).formula);
    triggers = {triggers{:}, addtriggers{:}};
end
% Check in reactions
for k=1:length(ms.reactions),
    addtriggers = gettrigger(ms.reactions(k).formula);
    triggers = {triggers{:}, addtriggers{:}};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD EVENTS WITH ABOVE TRIGGERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(triggers),
    if ~isempty(ms.states),
        sn = ms.states(1).name; % used for dummy assignments
    else
        sn = ms.parameters(1).name;
    end
    for k=1:length(triggers),
        ms.events(end+1).name = sprintf('piecewise_event_%d',k);
        ms.events(end).trigger = triggers{k};
        ms.events(end).assignment(1).variable = sn;
        ms.events(end).assignment(1).formula = sn;
        ms.events(end).notes = 'Just a dummy assignment for correct piecewise timing';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINISH IT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelout = SBmodel(ms);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helpfunction to get the trigger expressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trigger] = gettrigger(input)
trigger = {};
% find the start indices 
index = strfind(input,'piecewiseSB');
% get the triggers
for k=1:length(index),
    work = input;
    work = work(index(k)+12:end);
    popen = 1; offset = 1;
    while popen ~= 0,
        if work(offset) == '(',
            popen = popen + 1;
        elseif work(offset) == ')',
            popen = popen - 1;
        end
        offset = offset + 1;
    end
    work = work(1:offset-2);
    % explode elements
    terms = explodePCSB(work);
    trigger = {trigger{:} terms{2:2:end}};
end
return