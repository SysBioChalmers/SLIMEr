%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sol,model] = simulateGrowth(model,fluxData)
%
% Benjamín J. Sánchez. Last update: 2018-01-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol,model] = simulateGrowth(model,fluxData)

%Constrain all fluxes with exp. data
for i = 1:length(fluxData.rxnIDs)
    stdev = min([fluxData.stdevs(i),abs(fluxData.averages(i))/1.96]);
    LB    = fluxData.averages(i) - 1.96*stdev;      %C.I. of 95%
    UB    = fluxData.averages(i) + 1.96*stdev;      %C.I. of 95%
    model = changeRxnBounds(model,fluxData.rxnIDs(i),LB,'l');
    model = changeRxnBounds(model,fluxData.rxnIDs(i),UB,'u');
end

%Simulation 1: Maximize maintenance
posNGAM = strcmp(model.rxnNames,'non-growth associated maintenance reaction');
model   = changeRxnBounds(model,model.rxns(posNGAM),0,'l');
model   = changeRxnBounds(model,model.rxns(posNGAM),+1000,'u');
model   = changeObjective(model,model.rxns(posNGAM),+1);
sol     = optimizeCbModel(model);

%Simulation 2: Force NGAM, minimize sum(abs(fluxes))
obj   = sol.x(posNGAM);
model = changeRxnBounds(model,model.rxns(posNGAM),obj*0.999,'l');
model = changeRxnBounds(model,model.rxns(posNGAM),obj*1.001,'u');
sol   = optimizeCbModel(model,'min','one');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%