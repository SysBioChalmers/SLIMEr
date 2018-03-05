%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addLipidSpecies(model,metName,metFormula,exchange)
%
% Benjamín J. Sánchez. Last update: 2018-03-05
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
    model   = addReaction(model, ...                      %model
                          {['r_' newID],rxnName}, ...     %rxn
                          {metID}, ...                    %metabolites
                          -1, ...                         %stoichiometry
                          false, ...                      %reversibility
                          0, ...                          %LB
                          1000, ...                       %UB
                          0);                             %c
    printRxnFormula(model,['r_' newID],true,true,true);
end
                  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%