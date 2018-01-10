%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model,k] = scaleAbundancesInModel(model,data)
%
% Benjamín J. Sánchez. Last update: 2018-01-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,k] = scaleAbundancesInModel(model,data)

%Find optimal scaling factor:
save('temp.mat','model','data')
k0   = 1;
kOpt = fminsearch(@unusedLipid,k0);
save('temp.mat','model','data','kOpt')

%Find optimality range:
krange(1) = fminsearch(@(x) +minScaling(x),kOpt);
krange(2) = fminsearch(@(x) -minScaling(x),kOpt);
disp(['Optimality range: k = [ ' num2str(krange(1)) ' , ' num2str(krange(2)) ' ]'])

%Scale with the average of the range:
k     = mean(krange);
model = adjustModel(model,k,true);
disp(['Scaled lipid data in model: k = ' num2str(k)])
delete('temp.mat')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exchange = unusedLipid(k)

%Load model and adjust stoich coeffs of the tail pseudo-rxn:
temp  = load('temp.mat');
model = temp.model;
model = adjustModel(model,k,false);

%Optimize model:
data    = temp.data;
[sol,~] = simulateGrowth(model,data.fluxData);

%Objective function: unused tails or backbones
exchange_tails = sol.x(strcmp(model.rxnNames,'lipid - tails exchange'));
exchange_backs = sol.x(strcmp(model.rxnNames,'lipid - backbones exchange'));
exchange       = exchange_tails + exchange_backs;

disp(['Scaling abundance data: k = ' num2str(k) ' -> exchange = ' num2str(exchange)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k = minScaling(k)

%Load model and adjust stoich coeffs of the tail pseudo-rxn:
temp  = load('temp.mat');
model = temp.model;
model = adjustModel(model,k,true);

%Optimize model:
data    = temp.data;
[sol,~] = simulateGrowth(model,data.fluxData);

disp(['Finding scaling range: k = ' num2str(k) ' -> mu = ' num2str(sol.f)])

%Any unfeasible simulation returns the original value:
if isempty(sol.f)
    k = temp.kOpt;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%