%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% composition = getChainDistribution(model,chains)
%
% Benjamín J. Sánchez. Last update: 2017-12-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function composition = getChainDistribution(model,chains)

%Simulate model:
sol = optimizeCbModel(model);

%Find growth:
Xpos = strcmp(model.rxnNames,'growth');
mu   = sol.x(Xpos);

%Find chain distribution:
composition = zeros(length(chains),1);
for i = 1:length(chains)
    chainName = ['C' chains{i} ' chain [cytoplasm]'];
    chainPos  = strcmp(model.metNames,chainName);
    rxnPos    = model.S(chainPos,:) < 0;
    rxnStoich = abs(model.S(chainPos,rxnPos));
    
    %Update composition:
    composition(i) = sol.x(rxnPos)*rxnStoich/mu;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%