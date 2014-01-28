function [exp] = SBexperiment(varargin)
% SBexperiment: creates an experiment object defining experiment settings 
%
% USAGE:
% ======
% [exp] = SBexperiment()                creates an empty SBexperiment object
% [exp] = SBexperiment(SBstructure)     creates an SBexperiment object from a MATLAB
%                                       structure in the internal experiment format
% [exp] = SBexperiment(expin)           construction from a given SBexperiment object (expin)
% [exp] = SBexperiment('file.exp')      converting a experiment text description 
%                                       to an SBexperiment object.
%
% Output Arguments:
% =================
% exp: SBexperiment object 

% Information:
% ============
% Copyright (C) 2005-2007 Fraunhofer Chalmers Centre, Gothenburg, Sweden
% Main author: Henning Schmidt
% 
% Changes for the SBTOOLBOX2:
% 1/1/2008  Henning Schmidt, henning@sbtoolbox2.org
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


flag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK THE TYPE OF THE INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0,
    inputType = 'empty';
elseif nargin == 1 || nargin == 2,
    if strcmp('SBexperiment',class(varargin{1})),
        inputType = 'SBexperiment';
        expInput = varargin{1};
    elseif strcmp('struct',class(varargin{1})),
        inputType = 'SBstructure';
        SBstructure = varargin{1};
    elseif strcmp('char',class(varargin{1})),
        % check if '.txt' given as extension. If yes, then import text
        % description
        filename = varargin{1};
        if ~isempty(strfind(filename,'.exp')),
            inputType = 'TextExpFile';
        elseif strcmp('ExperimentAsTextString', varargin{2}),
            inputType = varargin{2};
        else
            error('Input argument of unknown type');
        end
    else 
        error('Input argument of unknown type');
    end
else
    error('Wrong number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT THE SBMODEL OBJECT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('empty',inputType),
    % Create empty SBstructure
    % parameter settings substructure
    paramicsettingsStruct = struct('name',{},'formula',{},'notes',{},'icflag',{});
    % parameter settings substructure
    parameterchangesStruct = struct('name',{},'formula',{},'notes',{});
    % event assignment substructure
    eventassignmentStruct = struct('variable',{},'formula',{});
    % state events substructure
    stateeventsStruct = struct('name',{},'trigger',{},'assignment',eventassignmentStruct,'notes',{});
    % Create SBstructure
    SBstructure = struct('name','unnamed_experiment','notes','no notes','paramicsettings',paramicsettingsStruct,'parameterchanges',parameterchangesStruct,'stateevents',stateeventsStruct);
    % construct the model object
    exp = class(SBstructure,'SBexperiment');
elseif strcmp('SBexperiment',inputType),
    % copy the model object
    exp = expInput;
elseif strcmp('SBstructure',inputType),
    % check if the given structure is a SBstructure (only check the
    % top-level fields)
    checkFields = {'name','notes','paramicsettings','parameterchanges','stateevents'};
    for k = 1:length(checkFields),
        if ~isfield(SBstructure,checkFields{k}),
            error('Given structure is not a valid internal SBexperiment structure.');
        end
    end
    % construct the model object
    exp = class(SBstructure,'SBexperiment');
elseif strcmp('TextExpFile',inputType),
    % check if a file with given filename exists
    [path,filename,ext,versn] = fileparts(filename);
    filename = fullfile(path, [filename '.exp']); 
    if ~exist(filename),
        errorMsg = sprintf('Experiment file, "%s", does not exist.', filename);
        error(errorMsg);
    end
    % If file exists then first load it
    expText = fileread(filename);
    % then convert it to SBstructure
    [SBstructure, errorMsg] = convertTextToExpSB(expText);
    % Check if error occurred while importing the SBML model
    if length(errorMsg) ~= 0,
        error(errorMsg);
    end
    % construct the model object
    exp = class(SBstructure,'SBexperiment');
elseif strcmp('ExperimentAsTextString', inputType),
    expText = varargin{1};
    % then convert text experiment to SBstructure
    [SBstructure, errorMsg] = convertTextToExpSB(expText);
    % Check if error occurred while importing the SBML model
    if length(errorMsg) ~= 0,
        error(errorMsg);
    end
    % construct the model object
    exp = class(SBstructure,'SBexperiment');
else
    error('Wrong input arguments.');
end
return
