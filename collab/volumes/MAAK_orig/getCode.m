function depcode = getCode(funName)

[list, builtins, classes, prob_files, prob_sym, eval_strings, called_from, java_classes] = depfun(funName,'-quiet');

nFiles=size(list,1);

depcode.list=list;
depcode.java_classes=java_classes;
depcode.builtins=builtins;
depcode.prob_files=prob_files;
depcode.called_from=called_from;

for i=1:nFiles
    callstr=['type(list{',num2str(i),'})'];
    depcode.mfiles(i).text=evalc(callstr);
end    