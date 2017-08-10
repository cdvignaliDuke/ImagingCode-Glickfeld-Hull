function normedTimeCourseEachTrial = getNormedTimeCourse(data,preEventFrames,...
    normalizeMethod)
% data is time-samples x trials
% normalizeMethod is either 'percent' or 'subtraction'
f0 = mean(data(1:preEventFrames,:),1);

if strcmp(normalizeMethod,'percent')
    normedTimeCourseEachTrial= bsxfun(@rdivide,data,f0);
elseif strcmp(normalizeMethod,'subtraction')
    normedTimeCourseEachTrial = bsxfun(@minus,data,f0);
end        
end
