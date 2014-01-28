function i = strmatchl(str,strarray);
%STRMATCHL Find loose matches for string. 
% I = STRMATCHL(STR,STRARRAY) Same as stringmatch but match can be
% anywhere within the string.
%
% 

if ~iscell(strarray)
    strarray = {strarray};
end

i = logical(zeros(size(strarray)));

for index = 1:length(strarray)
    
    if ~isempty(strfind(strarray{index},str))
        i(index)=true;
    end
end

return;