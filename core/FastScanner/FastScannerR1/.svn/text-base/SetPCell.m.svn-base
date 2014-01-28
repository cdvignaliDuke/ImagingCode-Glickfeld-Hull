function SetPCell
global fs

if fs.iniDone
    if isempty(fs.pCell.LUT)
        putsample(fs.DAQ.aoPockels,[fs.pCell.max fs.DAQ.aoCurrentTriggerValue]);
    else
       volt=Power2PockelsVoltage(fs.pCell.max);
       putsample(fs.DAQ.aoPockels,[volt fs.DAQ.aoCurrentTriggerValue]);
    end
end
