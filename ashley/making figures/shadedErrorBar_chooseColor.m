function figHandle = shadedErrorBar_chooseColor(x,y,yerr,lineColor)
%     use to set the patch color of shadedErrorBar to slightly lighter than
%     main line
    h = shadedErrorBar(x,y,yerr);
    h.mainLine.Color = lineColor;
    h.patch.FaceColor = errorPatchColor(h.mainLine.Color);
    h.edge(1).Color = 'none';
    h.edge(2).Color = 'none';
    figHandle = h;
end