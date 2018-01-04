%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sol = simulateGrowth(model)
%
% Benjamín J. Sánchez. Last update: 2017-12-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sol = simulateGrowth(model)

posG = strcmp(model.rxnNames,'D-glucose exchange');
posX = strcmp(model.rxnNames,'growth');

%Simulation 1: Limit glucose, maximize growth
model = changeRxnBounds(model,model.rxns(posG),-1,'b');
model = changeObjective(model,model.rxns(posX),+1);
sol   = optimizeCbModel(model);

%Simulation 2: Force growth, minimize sum(abs(fluxes))
mu    = sol.x(posX);
model = changeRxnBounds(model,model.rxns(posX),mu*0.999,'l');
sol   = optimizeCbModel(model,'min','one');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%