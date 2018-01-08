%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ATP,NADH,NADPH,netATP] = getEnergyexpenditure(model,flux)
%
% Benjamín J. Sánchez. Last update: 2018-01-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ATP,NADH,NADPH,netATP] = getEnergyexpenditure(model,flux)

%Get metabolites from all possible compartments:
ATPmets   = getMetFromAllComps(model,'ATP');
NADHmets  = getMetFromAllComps(model,'NADH');
NADPHmets = getMetFromAllComps(model,'NADPH');

%Add total energy consumption:
ATP   = 0;
NADH  = 0;
NADPH = 0;
for i = 1:length(model.rxns)
    netProduction = flux(i)*model.S(:,i);
    ATP   = ATP   + getConsumption(netProduction,ATPmets);
    NADH  = NADH  + getConsumption(netProduction,NADHmets);
    NADPH = NADPH + getConsumption(netProduction,NADPHmets);
end

%Get expenditure ratios:
ratioNADH  = getRatio(model,'NADH');
ratioNADPH = getRatio(model,'NADPH');

%Compute net ATP expenditure:
netATP = ATP + ratioNADH*NADH + ratioNADPH*NADPH;

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

function ratio = getRatio(model,met)

%Produce a 100% respiring model:
posNGAM = strcmp(model.rxnNames,'non-growth associated maintenance reaction');
sol_pre = simulateGrowth(model,model.rxns(posNGAM),'l');

%Force 1 mmol/gDWh of the met to go to waste:
if strcmp(met,'NADH')
    %           NADH   ->   NAD(+)  +   H+     (cytosolic)
    mets  = {'s_1203[c]','s_1198[c]','s_0794[c]'};
    model = addReaction(model,'NADHdrain',mets,[-1,+1,+1],false,+1,+1,0);
    %                                           stoich    rev   LB UB c
elseif strcmp(met,'NADPH')
    %           NADPH  ->  NADP(+)  +   H+     (cytosolic)
    mets  = {'s_1212[c]','s_1207[c]','s_0794[c]'};
    model = addReaction(model,'NADPHdrain',mets,[-1,+1,+1],false,+1,+1,0);
    %                                            stoich    rev   LB UB c
end

sol_post = simulateGrowth(model,model.rxns(posNGAM),'l');
ratio    = sol_pre.f - sol_post.f;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%