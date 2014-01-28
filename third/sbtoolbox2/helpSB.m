function [] = helpSB(varargin)
% helpSB: SBTOOLBOX2 help function displaying the web based user reference
% information using the MATLAB web browser. To be able to use this function
% you need to have an active internet connection.
%
% USAGE:
% ======
%           helpSB
%           helpSB 'SBT2 functionname'
%
% EXAMPLES:
% =========
%           helpSB
%           helpSB SBsimulate

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

if nargin == 0,
    web('http://www.sbtoolbox2.org/printview.php?display=documentationSBT&menu=userreference','-notoolbar');
elseif nargin == 1,
    if ischar(varargin{1}),
        web(sprintf('http://www.sbtoolbox2.org/printview.php?display=documentationSBT&menu=userreference#%s',varargin{1}),'-notoolbar');
    else
        web('http://www.sbtoolbox2.org/printview.php?display=documentationSBT&menu=userreference','-notoolbar');
    end
else
    web('http://www.sbtoolbox2.org/printview.php?display=documentationSBT&menu=userreference','-notoolbar');
end
