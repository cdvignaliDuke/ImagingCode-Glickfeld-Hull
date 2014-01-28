function out = tcRemoveDC(timeCourses)
% SCHEDULED FOR REMOVAL

    
[nframes,ncells] = size(timeCourses);

if any(isnan(timeCourses(:)))
	out = timeCourses-repmat(nanmean(timeCourses),nframes,1);
else
    out = timeCourses-repmat(mean(timeCourses),nframes,1);
end

return