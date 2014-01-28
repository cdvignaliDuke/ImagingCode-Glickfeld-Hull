function result = RegMcoreSlave(cmdstr, source, target)
%REGMCORESLAVE 

% global SRCTRACE;
% global RESULTTRACE;
% 
% SRCTRACE = [SRCTRACE squeeze(mean(mean(source,1),2))];
% clf
% plot(SRCTRACE);

sourceStack = ijarray2plus(source,'single'); 
targetStack = ijarray2plus(target,'single');
clear source;clear target;

al=IJAlign_AK;
disp(cmdstr);

resultStack = al.doAlign(cmdstr, sourceStack, targetStack);

clear sourceStack;clear targetStack;

% dummy object to disable automatic scaling in StackConverter
dummy = ij.process.ImageConverter(resultStack);
dummy.setDoScaling(0); % this static property is used by StackConverter

% stack = resultStack.getStack;
% ip = stack.getProcessor(1);
%disp(round([ip.getMin ip.getMax]));

converter = ij.process.StackConverter(resultStack); % don't use ImageConverter for stacks!
converter.convertToGray16;

% stack = resultStack.getStack;
% ip = stack.getProcessor(1);
%disp(round([ip.getMin ip.getMax]));

%source = ij2arrayAK(source);
result = ijplus2array(resultStack);

% this line moved to RegMultiCore
%source = permute(source,[1,3,2]); 

% RESULTTRACE = [RESULTTRACE squeeze(mean(mean(result,1),2))];
% 
% hold on;
% plot(RESULTTRACE,'r');
% 
% pause;
return;
