%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = scaleAbundancesInModel(model)
%
% Benjamín J. Sánchez. Last update: 2017-12-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = scaleAbundancesInModel(model)

%Find optimal scaling factor:
k0       = 1;
kOpt     = fminsearch(@unusedLipid,k0);
modelOpt = adjustModel(model,kOpt);

%Optimize model:
cd ../simulations
solOpt = simulateGrowth(modelOpt);
muOpt  = solOpt.f;
cd ../models
save('kOpt.mat','kOpt','muOpt')

%Find optimality range:
krange(1) = fminsearch(@(x) +minScaling(x),kOpt);
krange(2) = fminsearch(@(x) -minScaling(x),kOpt);
disp(['Optimality range: k = [ ' num2str(krange(1)) ' , ' num2str(krange(2)) ' ]'])

%Scale with the maximum of the range:
model = adjustModel(model,krange(2));
disp(['Scaled chain data in model: k = ' num2str(krange(2))])
delete('kOpt.mat')

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

function k = minScaling(k)

%Load model and adjust stoich coeffs of the tail pseudo-rxn:
model = load('yeast_7.8_SLIMEr.mat');
model = model.model_SLIMEr;
model = adjustModel(model,k);

%block exchange of tails and backbones:
posT = strcmp(model.rxnNames,'lipid - tails exchange');
posB = strcmp(model.rxnNames,'lipid - backbones exchange');
model = changeRxnBounds(model,model.rxns(posT),0,'b');
model = changeRxnBounds(model,model.rxns(posB),0,'b');

%Optimize model:
cd ../simulations
sol = simulateGrowth(model);
cd ../models

disp(['Finding scaling range: k = ' num2str(k) ' -> mu = ' num2str(sol.f)])

%Any simulation that results in impaired growth (>1%) returns the original value:
kOpt  = load('kOpt.mat');
muOpt = kOpt.muOpt;
if abs(sol.f - muOpt)/muOpt > 0.01
    k = kOpt.kOpt;
end

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