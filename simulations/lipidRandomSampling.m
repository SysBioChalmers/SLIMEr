%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% abundance = lipidRandomSampling(model,data,Nsim)
%
% Benjamín J. Sánchez. Last update: 2018-05-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function abundance = lipidRandomSampling(model,data,Nsim)

%Simulate a flux distribution with the corresponding model:
[sol,model] = simulateGrowth(model,data.fluxData);
posX          = strcmp(model.rxnNames,'growth');
muS           = sol.x(posX);                                %1/h

%Get a number of simulations from random sampling:
cd optGpSampler_1.1_Matlab
sModel  = optGpSampler(model,[],Nsim,500,4,'gurobi',1);
samples = sModel.points;
cd ..

%Find matching positions for each species and compute predicted abundance:
abundance = zeros(length(data.metNames),Nsim);
isSLIME   = contains(model.rxnNames,'SLIME rxn');
for j = 1:length(data.metNames)
    pos = matchToModel(model,data.metNames{j});
    if sum(pos) > 0
        %Get MW for each met:
        MWs = getMWfromFormula(model.metFormulas(pos));     %g/mmol
        
        %Get fluxes in which each met gets consumed for the SLIME rxn:
        isSub   = model.S(pos,:) < 0;
        fluxPos = isSub.*isSLIME';
        fluxes  = fluxPos*samples;
        
        %Some fluxes end up negative due to the random sampling algorithm:
        fluxes(fluxes < 0) = 0;                             %mmol/gDWh
        
        %Compute abundance predicted by model:
        abundance(j,:) = sum(fluxes.*MWs,1)/muS*1000;       %mg/gDW
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%