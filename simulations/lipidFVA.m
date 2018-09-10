%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [minVal,maxVal] = lipidFVA(model,lipidName,chainName)
%
% Benjamin J. Sanchez. Last update: 2018-09-03
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
    minVal = sum(model.c.*sol.x);
    
    %Maximize:
    sol    = optimizeCbModel(model,'max');
    maxVal = sum(model.c.*sol.x);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
