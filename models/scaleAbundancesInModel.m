%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model,k] = scaleAbundancesInModel(model)
%
% Benjamín J. Sánchez. Last update: 2018-01-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,k] = scaleAbundancesInModel(model)

%Find optimal scaling factor:
k0       = 1;
kOpt     = fminsearch(@unusedLipid,k0);
modelOpt = adjustModel(model,kOpt,true);

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

%Scale with the average of the range:
k     = mean(krange);
model = adjustModel(model,k,true);
disp(['Scaled lipid data in model: k = ' num2str(k)])
delete('kOpt.mat')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exchange = unusedLipid(k)

%Load model and adjust stoich coeffs of the tail pseudo-rxn:
model = load('yeast_7.8_SLIMEr.mat');
model = model.model_SLIMEr;
model = adjustModel(model,k,false);

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
model = adjustModel(model,k,true);

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