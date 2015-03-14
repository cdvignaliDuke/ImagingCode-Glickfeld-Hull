mouse = 'AW07';
SubNum = '607';
date = '150121';
ImgFolder = '003';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
%run FlahsingStim_dataStruct.m to organize data first

%%
% create structure dataMask for max dF/F for datasets: direction tuning,
% flashing stim (first flash rsp), flashing stim rsp to target (for each
% cycle)
% parameters for imCellEditInteractive are 0.9 and 3 pixels
siz = size(dataStructDFoverF.Cycles,2);

% direction tuning dataset
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\003'];
cd(CD);
load('mask&TCDir.mat');
dirTuning = mask_cell;
clear mask_cell

dataMasks.dirTuning = dirTuning;

% retinotopy dataset
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\001'];
cd(CD);
load('mask&TCRet.mat');
retTuning = mask_cell;
clear mask_cell

dataMasks.retTuning = retTuning;

% flashing stim datasets
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);

for icyc = 1:siz
    b = 5;
    data_reg = dataStructDFoverF.cycData{icyc};
    siz1 = size(data_reg);
    thisCorrMap = zeros(siz1(1),siz1(2));
    for ix = b:siz1(2)-b
        for iy = b:siz1(1)-b
            TC = data_reg(iy,ix,:);
            surround = (data_reg(iy-1,ix-1,:)+data_reg(iy-1,ix,:)+data_reg(iy-1,ix+1,:)+data_reg(iy,ix-1,:)+data_reg(iy,ix+1,:)+data_reg(iy+1,ix-1,:)+data_reg(iy+1,ix,:)+data_reg(iy+1,ix+1,:))/8;
            R = corrcoef(TC,surround);
            thisCorrMap(iy,ix) = R(1,2);
        end
    end
    corrMap{1,icyc} = thisCorrMap;
    clear data_reg siz1 surround R TC thisCorrMap
end


for icyc = 1:siz
    thisMax = corrMap{icyc};
    bwout = imCellEditInteractive(thisMax);
    thisMask = bwlabel(bwout);
    FScycCorrDF{icyc} = thisMask;
    clear thisMask thisMax bwout
end

dataMasks.FScycCorrDF = FScycCorrDF;

%flashing stim - average first response
%     start = 1;
% for icyc = 1:siz
%     thisData = dataStructDFoverF.cycData{icyc};
%     thisL = dataStructDFoverF.cycTrialL{icyc};
% 
%     for itrial = 1:dataStructDFoverF.cycNTrials{icyc}
%         thisDatabyTrial = thisData(:,:,(1+(thisL.*(itrial-1)):thisL.*itrial));
%         firstsData(:,:,start:start+29) = thisDatabyTrial(:,:,1:30);
%         start = start+30;
%     end
%     clear thisData thisL thisDatabyTrial 
% end
% avgFirstsData = mean(avgFirstsData,3);
% 
% 
% b = 5;
% data_reg = avgFirstsData;
% siz1 = size(data_reg);
% corrMapFirsts = zeros(siz1(1),siz1(2));
% for ix = b:siz1(2)-b
%     for iy = b:siz1(1)-b
%         TC = data_reg(iy,ix,:);
%         surround = (data_reg(iy-1,ix-1,:)+data_reg(iy-1,ix,:)+data_reg(iy-1,ix+1,:)+data_reg(iy,ix-1,:)+data_reg(iy,ix+1,:)+data_reg(iy+1,ix-1,:)+data_reg(iy+1,ix,:)+data_reg(iy+1,ix+1,:))/8;
%         R = corrcoef(TC,surround);
%         corrMapFirsts(iy,ix) = R(1,2);
%     end
% end
% clear data_reg surround R TC 
% 
% bwout = imCellEditInteractive(corrMapFirsts);
% FSfirstRsp = bwlabel(bwout);
% 
% dataMasks.FSfirstRsp = FSfirstRsp;
% 
% %flashing stim - target response
% for icyc = 1:siz
%     thisData = dataStructDFoverF.cycData{icyc};
%     thisL = dataStructDFoverF.cycTrialL{icyc};
%     start = 1;
%     for itrial = 1:dataStructDFoverF.cycNTrials{icyc}
%         thisDatabyTrial = thisData(:,:,(1+(thisL.*(itrial-1)):thisL.*itrial));
%         targetData(:,:,start:start+49) = thisDatabyTrial(:,:,thisL-49:thisL);
%         start = start+50;
%     end
%     cycTargetMax{icyc} = max(targetData,[],3);
%     clear thisData thisL thisDatabyTrial targetData
% end
% 
% for icyc = 1:siz
%     thisData = dataStructDFoverF.cycData{icyc};
%     thisL = dataStructDFoverF.cycTrialL{icyc};
%     start = 1;
%     for itrial = 1:dataStructDFoverF.cycNTrials{icyc}
%         thisDatabyTrial = thisData(:,:,(1+(thisL.*(itrial-1)):thisL.*itrial));
%         targetData(:,:,start:start+49) = thisDatabyTrial(:,:,thisL-49:thisL);
%         start = start+50;
%     end
%     clear thisData thisL thisDatabyTrial 
%     b = 5;
%     data_reg = targetData;
%     siz = size(data_reg);
%     thisCorrMap = zeros(siz(1),siz(2));
%     for ix = b:siz(2)-b
%         for iy = b:siz(1)-b
%             TC = data_reg(iy,ix,:);
%             surround = (data_reg(iy-1,ix-1,:)+data_reg(iy-1,ix,:)+data_reg(iy-1,ix+1,:)+data_reg(iy,ix-1,:)+data_reg(iy,ix+1,:)+data_reg(iy+1,ix-1,:)+data_reg(iy+1,ix,:)+data_reg(iy+1,ix+1,:))/8;
%             R = corrcoef(TC,surround);
%             thisCorrMap(iy,ix) = R(1,2);
%         end
%     end
%     corrMapTarget{1,icyc} = thisCorrMap;
%     clear data_reg siz surround R TC thisCorrMap targetData
% end
% 
% for icyc = 1:siz
%     thisMax = corrMapTarget{icyc};
%     bwout = imCellEditInteractive(thisMax);
%     thisMask = bwlabel(bwout);
%     FScycTargetRsp{icyc} = thisMask;
%     clear thisMask thisMax bwout
% end
% 
% dataMasks.FScycTargetRsp = FScycTargetRsp;

CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
save('dataMasks.mat','dataMasks')