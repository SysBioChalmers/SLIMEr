%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lipidRxns = getLipidRxns(model)
%
% Benjamín J. Sánchez. Last update: 2018-01-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ATP,NADH,NADPH] = getEnergyexpenditure(model,flux)

ATP   = 0;
NADH  = 0;
NADPH = 0;

ATPmets   = getMetFromAllComps(model,'ATP');
NADHmets  = getMetFromAllComps(model,'NADH');
NADPHmets = getMetFromAllComps(model,'NADPH');

for i = 1:length(model.rxns)
    netProduction = flux(i)*model.S(:,i);
    ATP   = ATP   + getConsumption(netProduction,ATPmets);
    NADH  = NADH  + getConsumption(netProduction,NADHmets);
    NADPH = NADPH + getConsumption(netProduction,NADPHmets);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mets = getMetFromAllComps(model,metName)

mets  = false(size(model.mets));
for i = 1:length(model.compNames)
    metName_i = [metName ' [' model.compNames{i} ']'];
    metPos_i  = strcmp(model.metNames,metName_i);
    if sum(metPos_i) > 0
        mets(metPos_i) = true;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function consumption = getConsumption(netProduction,mets)

metProduction = netProduction(mets);
consumption   = abs(metProduction(metProduction < 0));
consumption   = sum(consumption);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%