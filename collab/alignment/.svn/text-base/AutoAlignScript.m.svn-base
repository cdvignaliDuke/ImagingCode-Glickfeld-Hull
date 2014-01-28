%% File parameters
pn_vivo = 'G:\users\lindsey\analysisLG\100318\invivo_zoom1_RGB';
pn_vitro = 'G:\users\lindsey\analysisLG\100318\invitro_zoom5_RGB';
anal_pn = 'G:\users\lindsey\analysisLG\100318\alignment\';
zoom_vivo = 1;
zoom_vitro = 5;
size_vitro = 50;
resz_vivo = 2;
resz_vitro = 2;
pix_vivo = 256;
pix_vitro = 512;
resxy_vivo = (103*(zoom_vivo) + 15.25)./pix_vivo;
resxy_vitro = (691*(zoom_vitro^-1))./pix_vitro;

%% Find beads
invitro_B = findRedBids(pn_vitro, [resxy_vitro resxy_vitro resz_vitro]);
invivo_B = findRedBids(pn_vivo, [resxy_vivo resxy_vivo resz_vivo]);
save(fullfile(anal_pn, 'invitro_B.txt'), 'invitro_B','-ascii');
save(fullfile(anal_pn,'invivo_B.txt'), 'invivo_B', '-ascii');

%% Align beads and create initial stack
[transform_B rev_transform_B list_pairs_B dist_list_B]=run_D3align_L_2([anal_pn 'invivo_B.txt'], [anal_pn 'invitro_B.txt'],100000);
save([anal_pn 'coeffs_B'],'transform_B', 'rev_transform_B', 'list_pairs_B', 'dist_list_B');
[r_m_out_tr g_m_out_tr b_m_out_tr]=readTIFandtransform(pn_vivo, rev_transform_B,[pix_vitro pix_vitro size_vitro],[resxy_vivo resxy_vivo resz_vivo],[resxy_vitro resxy_vitro resz_vitro]);

%% Find cell pairs, align, and create final stack
save(fullfile(anal_pn, 'invitro_cells.txt'), 'invitro_cells','-ascii');
save(fullfile(anal_pn, 'invivo_transf_cells.txt'), 'invivo_transf_cells', '-ascii');
invivo_cells=transform(transform_B, invivo_transf_cells);
save(fullfile(anal_pn, 'invivo_cells.txt'), 'invivo_cells', '-ascii');
[transform rev_transform list_pairs dist_list]=run_D3align_L_2([anal_pn 'invivo_cells.txt'], [anal_pn 'invitro_cells.txt'],100000);
save([anal_pn 'coeffs'],'transform', 'rev_transform', 'list_pairs', 'dist_list');
[r_m_out_tr g_m_out_tr b_m_out_tr]=readTIFandtransform(pn_vivo, rev_transform,[pix_vitro pix_vitro size_vitro],[resxy_vivo resxy_vivo resz_vivo],[resxy_vitro resxy_vitro resz_vitro]);


