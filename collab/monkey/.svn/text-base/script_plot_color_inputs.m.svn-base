
% you need to load stats


contrast = [20.6; 8.4; 8.5; 8.8; 35.7; 9.7; 8.3; 7.3; 40; 80]/100;
Ncolors=10;

Ncells=length(stats);
resp=reshape([stats.dir_ratio_change], Ncolors, Ncells);
resp_normalized=resp./repmat(contrast,1,Ncells);
plot_color_inputs(resp_normalized);

filter=find([stats.p_value_resp]<0.05 );
plot_cone_contrast(resp,filter);
plot_cone_contrast(resp_normalized,filter);
