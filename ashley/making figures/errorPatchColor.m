function patchColor = errorPatchColor(lineColor)
%     use to set the patch color of shadedErrorBar to slightly lighter than
%     main line
    patchColor = lineColor +((1 - lineColor).*0.5);
end