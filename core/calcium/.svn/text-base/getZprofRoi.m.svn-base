function [xp yp] = getZprofRoi
%getZprofRoi (calcium): utility function to return roi from stackZProfiler window
%   getZProfRoi
%       prints x and y points in ROI polygon to command window
%   [xp,yp] = getZProfRoi
%       returns x and y points in output vectors
%
%   Extracts ROI information from current figure window
%
%$Id$

figH = gcf;

ud = get(figH, 'UserData');
if ~isfield(ud, 'xp')
    error('Current figure does not seem to be a stackZProfiler figure');
end

if nargout == 0
    % print it as a string
    fprintf(1, 'xp = %s;\n', mat2str(roundto(ud.xp(:)', 1)));
    fprintf(1, 'yp = %s;\n', mat2str(roundto(ud.yp(:)', 1)));
elseif nargout == 2
    xp = ud.xp; 
    yp = ud.yp;
else
    error('incorrect number of output arguments: should be 0 or 2');
end


    
