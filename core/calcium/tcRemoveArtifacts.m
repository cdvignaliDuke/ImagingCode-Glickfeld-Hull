function [result,coefs] = tcRemoveArtifacts(timeCourses,v)
%TCREMOVEARTIFACTS Removes artifacts from time courses
% [RESULT,COEFS]=TCREMOVEARTIFACTS(TIMECOURSES,ARTIFACTS).
% Projects ARTIFACTS onto TIMECOURSES one by one and subtract the result.
%
% by Vincent Bonin and Kenichi Ohki
%

% transform artifacts into orthogonal basis functions
% [nsamples,nartifacts]=size(artifacts);
% k = min(12,nartifacts);
% normalized = zscore(artifacts);
% [u,s,v]=svd(normalized','econ');

[nsamples,ncells]=size(timeCourses);
[temp,ncs]=size(v);

ind = 1:nsamples;

result = timeCourses ;

for ia = 1:ncs
    for ic = 1:ncells
        coefs(ia,ic) = (result(ind,ic)'*v(ind,ia)./norm(v(ind,ia))^2);
        result(:,ic)= result(:,ic)-coefs(ia,ic)*v(:,ia);
    end
end

ev = (var(timeCourses)-var(result))/var(timeCourses)*100;
%ev = sqrt(mean((result-timeCourses).^2))./var(timeCourses)*100;

fprintf('Artifacts explained %2.1f+/-%2.1f%% of the variance in time courses\n',mean(ev),std(ev));

return;
