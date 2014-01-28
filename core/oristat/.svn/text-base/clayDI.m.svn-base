function [best_dir_dp, null_dir_dp, R_best_dir_dp, R_null_dir_dp, DI_dp] = clayDI (dir_data, ori_vector_angle)

% Clay's direction index from dot product with average vector
% dir_data: 1 x Ndirections
% ori_vector_angle: preferred orientation angle
% CDI: Clay's direction index
% better_dir: better responding direction 
%
%   Kenichi Ohki 09/20/04
%

nstim_per_run = length(dir_data);

R1=0;
R2=0;

ave_vector = [cos(ori_vector_angle * pi /180); sin(ori_vector_angle * pi /180)];

for i=1:nstim_per_run
    dp= [cos(2*(i-1)*pi/nstim_per_run), sin(2*(i-1)*pi/nstim_per_run)] * ave_vector;
    if dp >= 0
        R1 = R1 + dp * dir_data(i);
    else
        R2 = R2 - dp * dir_data(i);
    end
end

if R1 >= R2
    best_dir_dp = ori_vector_angle;
    null_dir_dp = mod(ori_vector_angle+180, 360);
    R_best_dir_dp = R1;
    R_null_dir_dp = R2;
    DI_dp = 1-R2/R1;
else
    best_dir_dp = mod(ori_vector_angle+180, 360);
    null_dir_dp = ori_vector_angle;
    R_best_dir_dp = R2;
    R_null_dir_dp = R1;
    DI_dp = 1-R1/R2;
end

