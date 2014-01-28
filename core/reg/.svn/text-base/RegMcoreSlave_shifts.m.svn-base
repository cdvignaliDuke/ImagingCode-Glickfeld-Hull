function result = RegMcoreSlave_shifts(cmdstr, source, target)
%REGMCORESLAVE_shifts 

sourceStack = ijarray2plus(source,'single'); 
targetStack = ijarray2plus(target,'single');
clear source;clear target;

al=IJAlign_shifts;

result = al.doAlign(cmdstr, sourceStack, targetStack);

clear sourceStack;clear targetStack;