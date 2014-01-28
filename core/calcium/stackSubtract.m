function out = stackSubtract(array,frame)
%STACKSUBTRACT Subtract frame from stack array
%   OUT = STACKSUBTRACT(ARRAY,FRAME)

out = bsxfun(@minus,array,frame);

return