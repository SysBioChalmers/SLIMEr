%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = adjustModel(model,k,block)
%
% Benjamín J. Sánchez. Last update: 2017-12-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = adjustModel(model,k,block)

%Block exchange of tails and backbones:
if block
    posT  = strcmp(model.rxnNames,'lipid - tails exchange');
    posB  = strcmp(model.rxnNames,'lipid - backbones exchange');
    model = changeRxnBounds(model,model.rxns(posT),0,'b');
    model = changeRxnBounds(model,model.rxns(posB),0,'b');
end

%Find positions:
backRxn  = strcmp(model.rxnNames,'lipid pseudoreaction - backbone');
backMets = model.S(:,backRxn) < 0;

%Chain stoich coeffs:
model.S(backMets,backRxn) = k*model.S(backMets,backRxn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%