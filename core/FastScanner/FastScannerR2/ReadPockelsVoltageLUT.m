function PockelsVoltageLUT=ReadPockelsVoltageLUT(LUTfilename)
%function reads file which contains LUT for Pockels cell

%SY 03/02/04

% PockelsVoltageLUT - two column matrix, column 1- Pockels cell volts,
% column 2 - laser power [mW]


if nargin>0
    try
        PockelsVoltageLUTinit=load(LUTfilename,'-ascii');
        [cmax, indmax]=max(PockelsVoltageLUTinit(:,2));
        [cmin, indmin]=min(PockelsVoltageLUTinit(:,2));
        PockelsVoltageLUT=PockelsVoltageLUTinit(indmax:sign(indmin-indmax):indmin,1:2);
        PockelsVoltageLUT(:,2)=PockelsVoltageLUT(:,2)*1000; %convert from W to mW
    catch
        PockelsVoltageLUT=[];
        display(['Error loading ' LUTfilename]);
    end
else
        PockelsVoltageLUT=[];
end

