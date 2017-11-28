%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeChainComp(model,chainData)
%
% Benjamín J. Sánchez. Last update: 2017-11-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeChainComp(model,chainData)

%Figure out the IDs for all tails:
tailIDs = cell(size(chainData.metNames));
for i = 1:length(tailIDs)
    tailName   = [chainData.metNames{i} ' [cytoplasm]'];
    tailPos    = strcmp(model.metNames,tailName);
    tailIDs(i) = model.mets(tailPos);
end

%Define fields for new rxn:
newID   = getNewIndex(model.rxns);
rxnID   = ['r_' newID];
rxnName = 'lipid pseudoreaction - tail';
tailPos = strcmp(model.metNames,'lipid - tails [cytoplasm]');
tailID  = model.mets(tailPos);

%Add lipid pseudo-rxn for tails:
model = addReaction(model, ...                      %model
                    {rxnID,rxnName}, ...            %rxn
                    [tailIDs;tailID], ...           %metabolites
                    [-chainData.abundance;1], ...   %stoichiometry
                    false, ...                      %reversibility
                    0, ...                          %LB
                    1000, ...                       %UB
                    0);                             %c
printRxnFormula(model,rxnID,true,true,true);

%Add exchange rxn for tails:
newID   = getNewIndex(model.rxns);
rxnName = 'tail exchange';
model   = addReaction(model, ...                      %model
                      {['r_' newID],rxnName}, ...     %rxn
                      tailID, ...                     %metabolites
                      -1, ...                         %stoichiometry
                      false, ...                      %reversibility
                      0, ...                          %LB
                      1000, ...                       %UB
                      0);                             %c
printRxnFormula(model,['r_' newID],true,true,true);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%