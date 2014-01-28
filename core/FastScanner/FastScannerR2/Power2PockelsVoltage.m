function volt=Power2PockelsVoltage(power)
%function converts desired power of laser to
%voltage that should be applied to Pockels cell
%function uses calibration measurement stored in file

%SY 03/01/04
%mod for FastScanner - SY 03/27/07

global fs;
% fs.pCell.LUT - two column matrix, column 1- Pockels cell volts, column 2 - laser power

format compact;


if isempty(fs.pCell.LUT)
    display('No Pockels Cell LUT available');
    volt=power;
    if volt>2
        volt=2;
    end
    if volt<-2
        volt=-2;
    end
else
    if power<min(fs.pCell.LUT(:,2))
        power=min(fs.pCell.LUT(:,2))+0.000001;
        display(['Out of the power range, power set to min possible power: ' num2str(power)]);
    end
    if power>max(fs.pCell.LUT(:,2))
        power=max(fs.pCell.LUT(:,2));
        display(['Out of the power range, power set to max possible power : ' num2str(power)]);
    end
    try
    volt = interp1(fs.pCell.LUT(:,2),fs.pCell.LUT(:,1),power,'linear',-4);
    catch
        display('Error during LUT interpolation, attempt to use simple max-min interpolation');
        [cmax, indmax]=max(fs.pCell.LUT(:,2));
        [cmin, indmin]=min(fs.pCell.LUT(:,2));
    volt = interp1(fs.pCell.LUT([indmax indmin],2),fs.pCell.LUT([indmax indmin],1),power,'linear',-4);

    end
    if volt==-4
        display('Out of the power range, PockCell voltage set to 0');
        volt=0;
    end
end


