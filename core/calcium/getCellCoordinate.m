function [Centroid,EquivDiameter, Eccentricity,Extent,MajorAxisLength,MinorAxisLength] = getCellCoordinate(labelimg);

% return xy-coodinates of centers of cells (2 x Ncells)
% origin: top-left, x: left to right, y: top to bottom
%
% get also region statistics, usually for cell regions
% Uses regionprops in Image Processing toolbox
%  returns: Centroid (nLbl,2)
%    and 3 more arrays of nLbl elements ('Scalar'):
%  'EquivDiameter' - Scalar; the diameter of a circle with the same area as the 
%   region. Computed as sqrt(4*Area/pi). 
% 'Eccentricity' - Scalar; the eccentricity of the ellipse that has the same 
%   second-moments as the region. The eccentricity is the ratio of the distance 
%   between the foci of the ellipse and its major axis length. 
%   The value is between 0 and 1. (0 and 1 are degenerate cases; an ellipse 
%   whose eccentricity is 0 is actually a circle, while an ellipse whose eccentricity 
%   is 1 is a line segment.) 
%  'Extent' - Scalar; the proportion of the pixels in the bounding box that are also in the 
%    region.  Computed as the Area divided by area of the bounding box. CR:
%    this will be a small number for funny shapes
% 
% input: labelimage: 2D image with connected regions labelled with indices
% output:
% many variables with obvious names.
% CR started 040903

strr=regionprops(labelimg, 'Centroid', 'EquivDiameter', ...
        'Eccentricity','Extent','MajorAxisLength','MinorAxisLength');
Ncells=length(strr);
Centroid=reshape([strr.Centroid],2,Ncells);
%Centroid = cell2mat({strr.Centroid}')';
% this is a CR crazy kludge to make a 2 by nLbl array
    % by making an intermediate cell array....  otherwise I got a 1 by
    % 2*nLbl array
EquivDiameter = [strr.EquivDiameter]; 
Eccentricity = [strr.Eccentricity];
Extent = [strr.Extent ];
MajorAxisLength = [strr.MajorAxisLength];
MinorAxisLength  = [strr.MinorAxisLength];




