function [vector_angle, vector_mag, vector_tune] = vector_average(ydata)

% a set of parameters of direction (orientation) selectivity are obtained by vector averaging.
% vector_angle: preferred angle obtained by vector averaging, in degree (0-360).
% vector_mag: vector magnitude
% vector_tune: vector sum / scalar sum of singal changes to all directions,
%       according to Bonhoeffer et al. (1995)
% this function can be used to estimate parameters of both orientation & direction selectivity.
% when the signal change to one direction is nagative, it is replaced by zero,
% because vector_tune becomes strange. So, vector_tune parameter is between 0 and 1.
%
%
%   Kenichi Ohki 09/20/04
%


nstim_per_run = length(ydata);

Vx=0;
Vy=0;
sum=0;

for i=1:nstim_per_run
    a = max(ydata(i),0); 
    Vx=Vx+a*cos(2*(i-1)*pi/nstim_per_run);
    Vy=Vy+a*sin(2*(i-1)*pi/nstim_per_run);
    sum=sum+a;
end

vector_angle=mod(atan2(Vy,Vx)*180/pi,360);
vector_mag=(Vx^2+Vy^2)^(1/2);
if sum~=0; 
    vector_tune = vector_mag/sum;
else vector_tune = 0;
end