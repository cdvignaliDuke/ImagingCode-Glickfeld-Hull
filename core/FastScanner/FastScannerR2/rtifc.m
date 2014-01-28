function [x,map,details] = rtifc(filename,n)
%RTIFC Read a TIFF image.
%   [X,MAP,PHOTO] = READTIFC(FILENAME,N) reads the Nth image from the
%   TIFF file specified by the string variable FILENAME.  X contains the
%   image data; its class is uint8.  MAP contains the colormap
%   information if it exists; its class is uint16. MAP is empty if the
%   image is an RGB or grayscale; X is M-by-N-by-3 if the image is RGB.
%   PHOTO is the value of the PhotometricInterpretation TIFF tag for the
%   image.
%
%   [X,MAP,PHOTO] = READTIFC(FILENAME) reads the first image from the
%   file
%
%   See also WTIFC, IMWRITE, IMREAD, and IMFINFO.

%   Chris Griffin 6-12-96
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2005/11/15 01:13:17 $
%#mex

error('MATLAB:rtifc:missingMEX', 'Missing MEX-file RTIFC');

