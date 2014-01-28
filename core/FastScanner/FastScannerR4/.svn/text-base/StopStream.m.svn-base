function p=StopStream

global fs;

if isfield(fs,'StreamTimer') && ~isempty(fs.StreamTimer)
    stop(fs.StreamTimer);
    delete(fs.StreamTimer); fs.StreamTimer = [];
end

if ~fs.DAQ.DoNotSave && isfield(fs,'fid_curr') && ~isempty(fs.fid_curr)
    try 
        fclose(fs.fid_curr);
    catch
    end
end
if ~fs.DAQ.DoNotSave && isfield(fs,'fid_prev') && ~isempty(fs.fid_prev)
    try 
        fclose(fs.fid_prev);
    catch
    end
end
if isfield(fs,'map') && ~isempty(fs.map)
    try
        fclose(fs.map); 
    catch
    end
end

return;
