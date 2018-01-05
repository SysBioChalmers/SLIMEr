%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sol = simulateGrowth(model,objRxnName)
%
% Benjamín J. Sánchez. Last update: 2018-01-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sol = simulateGrowth(model,objRxnName)

posG = strcmp(model.rxnNames,'D-glucose exchange');
posR = strcmp(model.rxnNames,objRxnName);

%Simulation 1: Limit glucose, maximize growth
model = changeRxnBounds(model,model.rxns(posG),-1,'b');
model = changeObjective(model,model.rxns(posR),+1);
sol   = optimizeCbModel(model);

%Simulation 2: Force growth, minimize sum(abs(fluxes))
obj   = sol.x(posR);
model = changeRxnBounds(model,model.rxns(posR),obj*0.999,'l');
sol   = optimizeCbModel(model,'min','one');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%