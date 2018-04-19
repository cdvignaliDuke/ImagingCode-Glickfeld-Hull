function df_f_norm = normalizeTC(downsampled_movie)
[x,y,z2] = size(downsampled_movie);
movie2d = double(reshape(downsampled_movie, x*y, z2));

% normalize df/f
meanF = mean(movie2d, 2);

df_f = bsxfun(@minus, movie2d, meanF);

df_f = bsxfun(@rdivide, df_f, meanF);

df_f_norm = df_f;
% max_dff = max(df_f, [], 2);
% 
% df_f_norm = bsxfun(@rdivide, df_f, max_dff);
end