%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = changeLipidComp(model,lipidData)
%
% Benjamín J. Sánchez. Last update: 2018-03-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeLipidComp(model,lipidData)

%Detect if model has or not a met for backbones: If yes then create a new
%pseudo-rxn for backbones, if not then change the normal lipid pseudo-rxn.
metPos = strcmp(model.metNames,'lipid - backbones [cytoplasm]');
if sum(metPos) > 0
    newID   = getNewIndex(model.rxns);
    rxnID   = ['r_' newID];
    rxnName = 'lipid pseudoreaction - backbone';
    metID   = model.mets{metPos};
else
    rxnID   = 'r_2108';
    rxnName = 'lipid pseudoreaction';
    metID   = 's_1096[c]';
end

%Create lipid pseudo-rxn for backbones (or modify normal one):
model = addReaction(model, rxnID, ...
                    'reactionName', rxnName, ...
                    'metaboliteList', [lipidData.metIDs;metID], ...
                    'stoichCoeffList', [-lipidData.abundance;1], ...
                    'reversible', false, ...
                    'lowerBound', 0, ...
                    'upperBound', 1000);
printRxnFormula(model,rxnID,true,true,true);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%