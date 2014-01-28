function outStack = stack_rescale(inStack, inRange, outRange)
%$Id: stack_rescale.m 59 2007-12-11 01:00:07Z histed $

if nargin < 2 || isempty(inRange)
    inRange = [min(inStack(:)) max(inStack(:))];
end

if nargin < 3
    % use full output range
    outRange = [intmin(class(inStack)) intmax(class(inStack))];
end
 
% use integer arithmetic  (automatically if inStack is integer type)
scaleFactor = diff(outRange) ./ diff(inRange);
outStack = (inStack - inRange(1)) .* scaleFactor;


    
