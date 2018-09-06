%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sol,model] = simulateGrowth(model,fluxData)
%
% Benjamin J. Sanchez. Last update: 2018-09-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol,model] = simulateGrowth(model,fluxData)

%Constrain all fluxes with exp. data
for i = 1:length(fluxData.rxnIDs)
    %The variation should never be more than the abundance: 
    stdev = min([fluxData.stdevs(i),abs(fluxData.averages(i))/1.96]);
    %The variation should be at least a 5% (to avoid unfeasible simulations):
    stdev = max([stdev,abs(fluxData.averages(i))/1.96*0.05]);
    LB    = fluxData.averages(i) - 1.96*stdev;      %C.I. of 95%
    UB    = fluxData.averages(i) + 1.96*stdev;      %C.I. of 95%
    model = changeRxnBounds(model,fluxData.rxnIDs(i),LB,'l');
    model = changeRxnBounds(model,fluxData.rxnIDs(i),UB,'u');
end

%Simulation 1: Maximize maintenance
posMaint = strcmp(model.rxnNames,'non-growth associated maintenance reaction');
model    = changeRxnBounds(model,model.rxns(posMaint),0,'l');
model    = changeRxnBounds(model,model.rxns(posMaint),+1000,'u');
model    = changeObjective(model,model.rxns(posMaint),+1);
sol      = optimizeCbModel(model);

%Simulation 2: Force maintenance, minimize sum(abs(fluxes))
obj   = sol.x(posMaint);
model = changeRxnBounds(model,model.rxns(posMaint),obj*0.999,'l');
model = changeRxnBounds(model,model.rxns(posMaint),obj*1.001,'u');
sol   = optimizeCbModel(model,'min','one');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
