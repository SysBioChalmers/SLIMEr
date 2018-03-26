%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lipidValidation
%
% Benjamín J. Sánchez. Last update: 2018-03-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
addpath('../data')
addpath('../models')

for i = 1:length(model_corrComp_val)
    %Read experimental data for each condition:
    data = readEjsingData(i);
    data = convertEjsingData(data,model_SLIMEr_val{i},false);
    metNames      = data.metNames(1:end-1);         %filter out ergosterol
    abundance_exp = data.abundance(1:end-1)*1000;   %mg/gDW
    std           = data.std(1:end-1)*1000;         %mg/gDW
    
    %Simulate a flux distribution with the corresponding model:
    [sol,model] = simulateGrowth(model_SLIMEr_val{i},data.fluxData);
    posX        = strcmp(model.rxnNames,'growth');
    mu          = sol.x(posX);
    isSLIME     = contains(model.rxnNames,'SLIME rxn');
    
    %Find matching positions for each species and compute predicted abundance:
    abundance_mod = zeros(size(abundance_exp));
    for j = 1:length(metNames)
        pos = matchToModel(model,metNames{j});
        if sum(pos) > 0
            %Get MW for each met:
            MWs = getMWfromFormula(model.metFormulas(pos));     %g/mmol
            
            %Get fluxes in which each met gets consumed for the SLIME rxn:
            isSub   = model.S(pos,:) < 0;
            fluxPos = isSub.*isSLIME';
            fluxes  = fluxPos*sol.x;
            
            %Compute abundance predicted by model:
            abundance_mod(j) = sum(fluxes.*MWs)/mu*1000;	%mg/gDW
        end
    end
    
    %Plot data:
    barPlot(abundance_exp,metNames,'[mg/gDW]','r',25,1400,std);
    hold on
    plot(1:length(metNames),abundance_mod,'ob','markers',5)
    hold off
end
rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%