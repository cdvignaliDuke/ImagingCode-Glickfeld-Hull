function wtifc(varargin)
%WTIFC Write image data to a TIFF file.
%   WTIFC(RGB, [], FILENAME, COMPRESSION, DESCRIPTION)
%   writes the image in the uint8 array RGB to the
%   file specified by the string FILENAME, using the
%   compression option specified by the string COMPRESSION.
%   COMPRESSION may be 'none', 'rle', 'packbits', 'ccitt', 
%   'fax3', or 'fax4'.  The 'ccitt', 'fax3', and 'fax4' options
%   are for binary images only.  It defaults to 'none'.
%   DESCRIPTION may be any string; it defaults to the empty
%   string.
%
%   WTIFC(GRAY, [], FILENAME, COMPRESSION, DESCRIPTION)
%   writes a grayscale image.
%
%   WTIFC(X, MAP, FILENAME, COMPRESSION, DESCRIPTION)
%   writes an indexed image.
%
%   See also IMWRITE, IMREAD, and IMFINFO.

%   Copyright 1984-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2005/11/15 01:13:34 $
%#mex

error('MATLAB:wtifc:missingMEX', 'Missing MEX-file WTIFC');
