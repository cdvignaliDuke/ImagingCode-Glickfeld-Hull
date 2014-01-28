function [] = convert2SB(filename)
% convert2SB: This function converts the txt and txtbc format of the
% SBmodels from the first versions (SBTOOLBOX) to the new generation
% versions (SBTOOLBOX2).
%
% The only difference between the SBTOOLBOX and the SBTOOLBOX2 model 
% format is that in the TEXT and TEXTBC representation the additional 
% information (that is needed, e.g., for the SBML export) is enclosed
% in curled brackets, not in square brackets.
%
% This function will simply read in files ('filename'.txt or .txtbc) 
% in the old format, exchange square brackets against curly brackets 
% and save them as 'filename_2'.txt or txtbc.
% 
% Please be aware: The SBTOOLBOX2 model format allows the specification 
% of lookup tables. Here, square brackets are used. This means that 
% executing this function on SBTOOLBOX2 TEXT or TEXTBC files, containing
% lookup tables, is not advisable.

% split up filename
[path,filebasename,ext,versn] = fileparts(filename);
% do some checks
if (~strcmp(ext,'.txt') && ~strcmp(ext,'.txtc')),
    error('No valid file format (only .txt or .txtbc model files allowed).');
end
if ~isempty(path), 
    error('The model file to be converted needs to be in the current folder.');
end
if ~exist(filename),
    errorMsg = sprintf('Model file "%s" does not exist.', filename);
    error(errorMsg);
end
% read the file
modelText = fileread(filename);
% check if conversion necessary
if isempty(regexp(modelText,'[','once')) && isempty(regexp(modelText,']','once')),
    disp('No square brackets found => Not converting.')
    return
end
% exchange the parentheses
modelText = strrep(modelText,'[','{');
modelText = strrep(modelText,']','}');
% save the new model file
newfilename = strcat(filebasename,'_2',ext);
fid = fopen(newfilename,'w');
fwrite(fid,modelText);
fclose(fid);
% result message
disp('The model file has been converted to the new version.');


