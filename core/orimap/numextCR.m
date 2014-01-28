function numstr = numextCR(num,nchars)
numstr=num2str(floor(num));
nchartmp = length(numstr);
while nchartmp < nchars
	numstr=['0',numstr];
    nchartmp = length(numstr);
end
