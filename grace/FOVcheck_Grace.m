%% 
clear all
clear global
%% 
date = '200121';
mouse = 'i1317';
ImgFolder = '003';
run = '000';
doReg = 0;
nrun = size(ImgFolder,1);
zoom = '2.0';
mod = '200um';

CD = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\2P_Imaging\' mouse '\' date '_' mouse '\' ImgFolder];
cd(CD);
imgMatFile = [ImgFolder '_000_' run '.mat'];
load(imgMatFile);

nframes = info.config.frames;
data = sbxread([ImgFolder '_000_' run],0,nframes);
fprintf(['Loaded ' num2str(nframes) ' frames \r\n'])

%Register Data, use Green PMT
if size(data,1) > 1
    data = data(1,:,:,:);
end
data = squeeze(data);

data_avg = mean(data(:,:,400:450),3);
[out data_reg] = stackRegister(data,data_avg);
data_reg_avg = mean(data_reg,3);


%Making the Figure
figure;  imagesc(data_reg_avg);
axis off
% title([date ' zoom' zoom ' ' mod ' ' mouse ' 920 Green Channel'])
truesize

% if exist(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_FOV_Check']);
%     print(['Z:\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\'  date '_' mouse '_FOV_Check\' date '_' mouse '_run' run '_zoom' zoom '_' mod '_FOVgreen.pdf'], '-dpdf','-bestfit')
%     save(['Z:\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_FOV_Check\' date '_' mouse '_run' run '_zoom' zoom '_' mod '_avgFOVgreen.mat'],'data_reg_avg')
%     writetiff(data_reg_avg, ['Z:\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_FOV_Check\' date '_' mouse '_run' run '_zoom' zoom '_' mod '_avgFOVgreen.tiff'])
%     
% else mkdir(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_FOV_Check']);
%     print(['Z:\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\'  date '_' mouse '_FOV_Check\' date '_' mouse '_run' run '_zoom' zoom '_' mod '_FOVgreen.pdf'], '-dpdf','-bestfit')
%     save(['Z:\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_FOV_Check\' date '_' mouse '_run' run '_zoom' zoom '_' mod '_avgFOVgreen.mat'],'data_reg_avg')
%     writetiff(data_reg_avg, ['Z:\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_FOV_Check\' date '_' mouse '_run' run '_zoom' zoom '_' mod '_avgFOVgreen.tiff'])
% end
% 
% 
