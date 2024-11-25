function [bnds,model] = PrepareModel(model)
% This function creates some necessary fields for the master file,
% like model.qualDir which indicates if a reaction is bidirectional or not.
[selExc, ~] = findExcRxns(model);
exchange=model.rxns(selExc);
selExcIds=find(selExc);
exclb=model.lb(selExcIds);
Exclb=num2cell(exclb);
excub=model.ub(selExcIds);
Excub=num2cell(excub);
bnds=[exchange Exclb Excub];

model.rxnID = findRxnIDs(model,model.rxns);
model.metIDs = findMetIDs(model,model.mets);

 minmax = [model.lb,model.ub];
 
model.qualDir(minmax(:,1)>=0 & minmax(:,2)>0)=1;
model.qualDir(minmax(:,1)<0 & minmax(:,2)<=0)=-1;
model.qualDir(minmax(:,1)<0 & minmax(:,2)>0)=0;
end 