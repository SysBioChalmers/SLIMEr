%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = adjustModel(model,k)
%
% Benjam�n J. S�nchez. Last update: 2017-12-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = adjustModel(model,k)

%Find positions:
backRxn  = strcmp(model.rxnNames,'lipid pseudoreaction - backbone');
backMets = model.S(:,backRxn) < 0;

%Chain stoich coeffs:
model.S(backMets,backRxn) = k*model.S(backMets,backRxn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%