function rng = scatterAxisRange(xData,yData)

minRng = min(min(xData),min(yData));
maxRng = max(max(xData),max(yData));
rng = [minRng maxRng];

end