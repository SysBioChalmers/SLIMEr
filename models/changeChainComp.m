%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeChainComp(model,chainData)
%
% Benjamín J. Sánchez. Last update: 2018-05-20
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
model = addReaction(model, rxnID, ...
                    'reactionName', rxnName, ...
                    'metaboliteList', [tailIDs;tailID], ...
                    'stoichCoeffList', [-chainData.abundance;1], ...
                    'reversible', false, ...
                    'lowerBound', 0, ...
                    'upperBound', 1000);
                
printRxnFormula(model,rxnID,true,true,true);
model.rxnConfidenceScores(strcmp(model.rxns,rxnID)) = 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%