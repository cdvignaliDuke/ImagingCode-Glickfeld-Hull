%find the mean RT for a given mouse across days. 
%find when tft was removed
clear
days = {'150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41'}; %'150718_img27', '150719_img27',
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
BEHAVE_DATABASE = '\\crash\data\home\andrew\Behavior\Data\';
BX_OUTPUT_DIR = 'Z:\Analysis\LeverAnalysis\BxAndAnalysisOutputs\BxOutputs\';
for kk=1:length(days)
   day_name  =  days{kk};
   bfile = dir([BEHAVE_DIR 'data-*i9' days{kk}(end-1:end) '-' days{kk}(1:6) '*' ]);
   behave_dest = [BEHAVE_DIR bfile.name];
   assert(length(bfile)) =1;
   b_data = load(behave_dest);
   TFTMat=[];
   bx_date = str2num(days{kk}(1:6));
   for i = 1:200;
       curr_bfile = dir([BEHAVE_DATABASE 'data-*i9' days{kk}(end-1:end) '-' num2str(bx_date-i) '*' ]);
       if exist([BEHAVE_DATABASE curr_bfile.name])==2;
           sub_bfile = dir([BEHAVE_DATABASE 'data-*i9' days{kk}(end-1:end) '-' num2str(bx_date-i) '*' ]);
           curr_bx = load([BEHAVE_DATABASE sub_bfile.name]);
           currTFT =  curr_bx.input.tooFastTimeMs;
           reactTimesMs = double(cell2mat(curr_bx.input.reactTimesMs));
           reactTimesMs(reactTimesMs<150)=[];
           meanCorrRT = mean(reactTimesMs);
           stdCorrRT = std(reactTimesMs);
           TFTMat = [TFTMat, currTFT];
       end
       if length(TFTMat)==10;
           break
       end
   end
   figure; plot(TFTMat)
   title([days{kk}, 'Too Fast Time setting for previous 10 days']);
   ylabel('time (ms)');
   xlabel('days');
end

