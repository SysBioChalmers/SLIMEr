%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addLipidSpecies(model,metName,metFormula)
%
% Benjam�n J. S�nchez. Last update: 2017-11-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addLipidSpecies(model,metName,metFormula)

newID = getNewIndex(model.mets);
metID = ['s_' newID '[c]'];
model = addMetabolite(model,metID, ...
                      'metName',[metName ' [cytoplasm]'], ...
                      'metFormula',metFormula);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%