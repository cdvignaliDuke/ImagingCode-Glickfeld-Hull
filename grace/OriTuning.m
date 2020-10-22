%% get path names
date = '200114';
ImgFolder = strvcat('003');
time = strvcat('1217');
mouse = 'i1314';
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
gl_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\jerry\2P_Imaging_Grace';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P';
load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningAndFits.mat']));
load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']));
load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']));

nCells = size(npSub_tc,2);
dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
ori_mat = dir_mat;
ori_mat(find(dir_mat>=180)) = ori_mat(dir_mat>=180)-180;
oris = unique(ori_mat);
for ic = 20
    figure;
    errorbar(oris,avgResponseEaOri(ic,:), semResponseEaOri(ic,:),'k-o','LineStyle','none','MarkerEdgeColor','k')
    hold on
    plot(0:180,vonMisesFitAllCellsAllBoots(:,1,ic),'b');
    figXAxis([ ],'Orientation (Deg)',[0 180],0:90:180,0:90:180);
    ylabel('dF/F','FontSize');
%     title(ic);
end
mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['OriTuning']));
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['OriTuning'], ['_oriTuning_' ic '.pdf']),'-dpdf')
