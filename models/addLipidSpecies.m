%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addLipidSpecies(model,metName,metFormula,exchange)
%
% Benjamín J. Sánchez. Last update: 2018-03-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addLipidSpecies(model,metName,metFormula,exchange)

newID    = getNewIndex(model.mets);
compName = metName(strfind(metName,'[')+1:strfind(metName,']')-1);
compPos  = strcmp(model.compNames,compName);
compID   = model.comps{compPos};
metID    = ['s_' newID '[' compID ']'];
model    = addMetabolite(model,metID, ...
                         'metName',metName, ...
                         'metFormula',metFormula);

%Add exchange rxn if indicated:
if exchange
    newID   = getNewIndex(model.rxns);
    metName = metName(1:strfind(metName,'[')-2);
    rxnName = [metName ' exchange'];
    model   = addReaction(model,['r_' newID], ...
                          'reactionName', rxnName, ...
                          'metaboliteList', {metID}, ...
                          'stoichCoeffList', -1, ...
                          'reversible', false, ...
                          'lowerBound', 0, ...
                          'upperBound', 1000);
                      
    printRxnFormula(model,['r_' newID],true,true,true);
end
                  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%