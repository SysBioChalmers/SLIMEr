%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = adjustModel(model,k,block,scaling)
%
% Benjamín J. Sánchez. Last update: 2018-03-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = adjustModel(model,k,block,scaling)

%Block exchange of tails and backbones:
if block
    posT  = strcmp(model.rxnNames,'lipid - tails exchange');
    posB  = strcmp(model.rxnNames,'lipid - backbones exchange');
    model = changeRxnBounds(model,model.rxns(posT),0,'b');
    model = changeRxnBounds(model,model.rxns(posB),0,'b');
end

%Switch what to rescale depending on flag:
switch scaling
    case 'backbones'
        rxnName = 'lipid pseudoreaction - backbone';
    case 'tails'
        rxnName = 'lipid pseudoreaction - tail';
end

%Find positions:
scaleRxn  = strcmp(model.rxnNames,rxnName);
scaleMets = model.S(:,scaleRxn) < 0;

%Change stoich coeffs:
model.S(scaleMets,scaleRxn) = k*model.S(scaleMets,scaleRxn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%