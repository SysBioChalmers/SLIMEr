%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sol = simulateGrowth(model,objRxn,glucOpt)
%
% Benjamín J. Sánchez. Last update: 2018-01-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sol = simulateGrowth(model,objRxn,glucOpt)

posG = strcmp(model.rxns,'1714');       %D-glucose exchange
posR = strcmp(model.rxns,objRxn);       %Objective rxn

%Simulation 1: Limit glucose, maximize growth
model = changeRxnBounds(model,model.rxns(posG),-1,glucOpt);
model = changeRxnBounds(model,model.rxns(posR),0,'l');
model = changeRxnBounds(model,model.rxns(posR),+1000,'u');
model = changeObjective(model,model.rxns(posR),+1);
sol   = optimizeCbModel(model);

%Simulation 2: Force growth, minimize sum(abs(fluxes))
obj   = sol.x(posR);
model = changeRxnBounds(model,model.rxns(posR),obj*0.999,'l');
sol   = optimizeCbModel(model,'min','one');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%