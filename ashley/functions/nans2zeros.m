function arrayWithZeros = nan2zeros(arrayWithNans)
arraySize = size(arrayWithNans);
arrayWithZeros = arrayWithNans;
arrayWithZeros(isnan(arrayWithNans)) = 0;
end