function [beta, t_value, variance_matrix]=regression_analysis (t, s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spm-like multiple regression analysis
%
%   2009. 5. 13.        K. Ohki
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    s.Non=10;
    s.Noff=10;
    s.Nstim_per_run=12;
    s.Nrep=8;
end

tau=1;
frame_in_sec=0.4;
lowcut_frame = 120;
%%%%%%%%%%%%%%%%%%%%%%%%

nframes_per_run = (s.Non+s.Noff)*s.Nstim_per_run;

Nframes = nframes_per_run * s.Nrep;

norm_t=t./mean(t)-1;
Ca_ker=Calcium_kernel(tau,frame_in_sec);

% smooth timecourse by Ca kernel
K=Toeplitz_matrix(Ca_ker,Nframes);
Kt=K*norm_t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiple regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% regressors
regressors = make_regressors (s, lowcut_frame, Ca_ker);
Nregressors=size(regressors,2);
sm_regressors = K * regressors;

figure;
imagesc(sm_regressors);
colormap(gray);
title('regressors');

% regression

beta=sm_regressors\Kt;
figure;
plot(beta(1:s.Nstim_per_run));
title('regression coefficients (beta)');

estimated_t = regressors*beta;
estimated_resp = regressors(:,1:s.Nstim_per_run)*beta(1:s.Nstim_per_run);
estimated_slow = regressors(:,s.Nstim_per_run+1:Nregressors)*beta(s.Nstim_per_run+1:Nregressors);

estimated_Kt = sm_regressors*beta;

figure;
subplot(2,1,1);
plot([norm_t, estimated_t]);
title('estimated timecourses');
subplot(2,1,2);
plot([estimated_resp, estimated_slow]);
title('visual vs slow modulation');

t_avg=mean(reshape(norm_t,nframes_per_run,s.Nrep),2);
estimated_t_avg=mean(reshape(estimated_t,nframes_per_run,s.Nrep),2);
estimated_resp_avg=mean(reshape(estimated_resp,nframes_per_run,s.Nrep),2);
estimated_slow_avg=mean(reshape(estimated_slow,nframes_per_run,s.Nrep),2);
figure
subplot(2,1,1);
plot([t_avg, estimated_t_avg]);
title('average timecourses');
subplot(2,1,2);
plot([estimated_resp_avg, estimated_slow_avg]);
title('average visual vs slow modulation');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error_t=Kt-estimated_Kt;
error_sum = error_t'*error_t;

V=K*K';
R = K - sm_regressors*pinv(sm_regressors)*K;

sigma = (error_sum / trace(R*V))^(1/2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t-value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KX=sm_regressors;
variance_matrix = inv(KX'*KX)*KX'*V*KX*inv(KX'*KX) .*(sigma^2);

t_value = zeros(s.Nstim_per_run,1);
for i=1:s.Nstim_per_run
    t_value(i) = beta(i)/(variance_matrix(i,i)^(1/2));
end

figure;
plot(t_value);
title('t values');
    


