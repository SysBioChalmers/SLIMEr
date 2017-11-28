%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addLipidSpecies(model,metName)
%
% Benjamín J. Sánchez. Last update: 2017-11-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addLipidSpecies(model,metName)

newID = getNewIndex(model.mets);
metID = ['s_' newID '[c]'];
model = addMetabolite(model,metID,'metName',[metName ' [cytoplasm]']);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%