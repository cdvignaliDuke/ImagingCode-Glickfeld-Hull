function imageStackNorm = normalizeImageStack(imageStack)
    imageStackZeroed = imageStack - min(imageStack(:));
    normPixelValue = max(imageStack(:));
    imageStackNorm = imageStackZeroed./normPixelValue;
end