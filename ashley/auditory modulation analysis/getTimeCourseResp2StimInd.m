function timeCourseResp2StimIndEaCond = getTimeCourseResp2StimInd(...
    neuronTimeCourseData,timeCourseField,stimInd,varargin)

if isempty(varargin)
    conditions = 1:size(neuronTimeCourseData,2);
else
    conditions = varargin{1};
end
nconditions = length(conditions);

timeCourseResp2StimIndEaCond = cell(1,nconditions);
for icond = conditions
    thisConditionData = eval(['neuronTimeCourseData(icond).' timeCourseField]);
    [timeCourseLength,ncells,~] = size(thisConditionData);
    
    timeCourseAtInd = nan(timeCourseLength,ncells);
    for icell = 1:ncells
        thisCellInd = stimInd(icell);
        thisCellTimeCourse = thisConditionData(:,icell,thisCellInd);
        timeCourseAtInd(:,icell) = thisCellTimeCourse;
    end
timeCourseResp2StimIndEaCond{icond} = timeCourseAtInd;
end

end