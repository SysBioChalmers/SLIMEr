%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addLipidSpecies(model,metName,metFormula,exchange)
%
% Benjamín J. Sánchez. Last update: 2017-12-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addLipidSpecies(model,metName,metFormula,exchange)

newID = getNewIndex(model.mets);
metID = ['s_' newID '[c]'];
model = addMetabolite(model,metID, ...
                      'metName',[metName ' [cytoplasm]'], ...
                      'metFormula',metFormula);

%Add exchange rxn if indicated:
if exchange
    newID   = getNewIndex(model.rxns);
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