%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [minVal,maxVal] = lipidFVA(model,lipidName,chainName)
%
% Benjamín J. Sánchez. Last update: 2017-12-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [minVal,maxVal] = lipidFVA(model,lipidName,chainName)

%Find metabolite positions:
lipidPos = ~cellfun(@isempty,strfind(model.metNames,[lipidName ' [']));
chainPos = strcmp(model.metNames,chainName);

%Find rxns in which the lipid species and the chain are formed:
lipidRxns = sum(model.S(lipidPos,:) > 0,1);
chainRxns = model.S(chainPos,:) > 0;
rxnPos    = lipidRxns.*chainRxns == 1;
rxnIDs    = model.rxns(rxnPos);

%Assign coefficients to objective function:
objCoeffs = model.S(chainPos,rxnPos)';

%Compute variability only when rxns exist:
if sum(objCoeffs) == 0
    minVal = 0;
    maxVal = 0;
else
    %Change objective function:
    model = changeObjective(model,rxnIDs,objCoeffs);

    %Minimize:
    sol    = optimizeCbModel(model,'min');
    minVal = sol.f;
    
    %Maximize:
    sol    = optimizeCbModel(model,'max');
    maxVal = sol.f;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%