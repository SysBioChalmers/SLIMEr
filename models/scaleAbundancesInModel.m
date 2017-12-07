%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = scaleAbundancesInModel(model)
%
% Benjamín J. Sánchez. Last update: 2017-12-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = scaleAbundancesInModel(model)

k0 = 1;
k  = fminsearch(@unusedLipid,k0);
disp(['Scaled chain data in model: k = ' num2str(k)])
model = adjustModel(model,k);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exchange = unusedLipid(k)

%Load model and adjust stoich coeffs of the tail pseudo-rxn:
model = load('yeast_7.8_SLIMEr.mat');
model = model.model_SLIMEr;
model = adjustModel(model,k);

%Optimize model:
cd ../simulations
sol = simulateGrowth(model);
cd ../models

%Objective function: unused tails or backbones
exchange_tails = sol.x(strcmp(model.rxnNames,'lipid - tails exchange'));
exchange_backs = sol.x(strcmp(model.rxnNames,'lipid - backbones exchange'));
exchange       = exchange_tails + exchange_backs;

disp(['Scaling abundance data: k = ' num2str(k) ' -> exchange = ' num2str(exchange)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = adjustModel(model,k)

%Find positions:
chainRxn  = strcmp(model.rxnNames,'lipid pseudoreaction - tail');
chainMets = model.S(:,chainRxn) < 0;

%Chain stoich coeffs:
model.S(chainMets,chainRxn) = k*model.S(chainMets,chainRxn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%