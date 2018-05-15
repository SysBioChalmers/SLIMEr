%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [lipids,cost] = compareCosts(model)
%
% Benjamín J. Sánchez. Last update: 2018-05-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lipids,cost] = compareCosts(model)

%Provide glucose & block growth:
posGluc   = strcmp(model.rxnNames,'D-glucose exchange');
posGrowth = strcmp(model.rxnNames,'growth');
model     = changeRxnBounds(model,model.rxns(posGluc),-1,'l');
model     = changeRxnBounds(model,model.rxns(posGrowth),0,'b');

%Cycle through all generic species in model:
LPRpos = strcmp(model.rxnNames,'lipid pseudoreaction - backbone');
lipPos = find(model.S(:,LPRpos) < 0);
lipids = cell(200,1);
cost   = zeros(200,2);
k = 1;
for i = 1:length(lipPos)
    if ~contains(model.metNames{lipPos(i)},'ergosterol [')
        %Find specific positions and go one level up if they are not SLIME rxns:
        specPos = find(model.S(lipPos(i),:) > 0);
        if ~contains(model.rxnNames{specPos(1)},'SLIME')
            lipPos(i) = find(model.S(:,specPos(1)) < 0);
            specPos   = find(model.S(lipPos(i),:) > 0);
        end
        
        %Allow free exchange of the generic species:
        model_i  = addReaction(model,'genOut','metaboliteList',model.mets(lipPos(i)), ...
                               'stoichCoeffList',-1,'lowerBound',0,'upperBound',1000);
        rxnPos   = strcmp(model_i.rxns,'genOut');
        
        %Cycle through all specific positions:
        for j = 1:length(specPos)
            %Analyze growth:
            rxnName   = model_i.rxnNames{specPos(j)};
            lipids{k} = rxnName(1:strfind(rxnName,' [')-1);
            [cost(k,1),cost(k,2)] = analyzeGRowth(model_i,specPos(j),rxnPos);
            k = k + 1;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [molCost,massCost] = analyzeGRowth(model,specPos,genPos)

%Optimize for the given rxn:
model    = changeObjective(model,model.rxns(specPos),+1);
sol      = optimizeCbModel(model);
molCost  = sol.f;
massCost = sol.x(genPos);
disp([model.rxnNames{specPos} ': ' num2str(molCost) ' mmol/gDWh = ' num2str(massCost) ' g/gDWh'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%