function h = scatterAVResponses(visResponse,audResponse,responseLim)
h = scatter(visResponse,audResponse,'ko');
h.MarkerFaceColor = 'k';
h.MarkerEdgeColor = [1,1,1];
hold on
plot(responseLim,responseLim,'k--');
figXAxis([],'visual',responseLim);
figYAxis([],'auditory',responseLim)
figAxForm([])
end