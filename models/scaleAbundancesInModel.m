%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model,k] = scaleAbundancesInModel(model,data,scaling)
%
% Benjamin J. Sanchez. Last update: 2018-09-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,k] = scaleAbundancesInModel(model,data,scaling)

%Find optimal scaling factor:
k0   = 1;
kOpt = fminsearch(@(k)unusedLipid(k,model,data,scaling),k0);

%Find optimality range:
krange(1) = fminsearch(@(k) +minScaling(k,model,data,scaling,kOpt),kOpt);
krange(2) = fminsearch(@(k) -minScaling(k,model,data,scaling,kOpt),kOpt);
disp(['Optimality range: k = [ ' num2str(krange(1)) ' , ' num2str(krange(2)) ' ]'])
if(krange(1) == krange(2))
    error('Could not find an optimality range!')
end

%Scale with the average of the range:
k     = mean(krange);
model = adjustModel(model,k,true,scaling);
disp(['Scaled ' scaling(1:end-1) ' data in model: k = ' num2str(k)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exchange = unusedLipid(k,model,data,scaling)

%Adjust stoich coeffs of the corresponding pseudo-rxn:
model = adjustModel(model,k,false,scaling);

%Optimize model:
try
    [sol,~] = simulateGrowth(model,data.fluxData);
catch
    sol.x = ones(size(model.rxns));
end

%Objective function: unused tails or backbones
exchange_tails = sol.x(strcmp(model.rxnNames,'lipid - tails exchange'));
exchange_backs = sol.x(strcmp(model.rxnNames,'lipid - backbones exchange'));
exchange       = exchange_tails + exchange_backs;

disp(['Scaling abundance data: k = ' num2str(k) ' -> exchange = ' num2str(exchange)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k = minScaling(k,model,data,scaling,kOpt)

%Adjust stoich coeffs of the corresponding pseudo-rxn:
model = adjustModel(model,k,true,scaling);

%Optimize model:
try
    [sol,~] = simulateGrowth(model,data.fluxData);
    posNGAM = strcmp(model.rxnNames,'non-growth associated maintenance reaction');
    disp(['Finding scaling range: k = ' num2str(k) ' -> Maintenance = ' num2str(sol.x(posNGAM))])
catch
    disp(['Finding scaling range: k = ' num2str(k) ' -> Maintenance = ' num2str(0)])
    k = kOpt;  %any unfeasible simulation returns the original value
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
