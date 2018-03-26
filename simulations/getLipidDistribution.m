%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = getLipidDistribution(model,lipidNames,chains,fluxData)
%
% Benjamín J. Sánchez. Last update: 2018-03-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = getLipidDistribution(model,lipidNames,chains,fluxData)

%Simulate model:
[sol,model_cons] = simulateGrowth(model,fluxData);

%Find growth:
posGluc  = strcmp(model.rxnNames,'D-glucose exchange');
posX     = strcmp(model.rxnNames,'growth');
posMaint = strcmp(model.rxnNames,'non-growth associated maintenance reaction');
mu       = sol.x(posX);

for i = 1:length(chains)
    chains{i} = ['C' chains{i} ' chain [cytoplasm]'];
end

composition = zeros(length(lipidNames),length(chains));

%Go through all SLIME rxns to find abundances:
SLIMEpos = find(contains(model.rxnNames,'SLIME rxn'));
for i = 1:length(SLIMEpos)
    %Find flux and all metabolites produced in each SLIME rxn:
    flux      = sol.x(SLIMEpos(i));
    metPos    = model.S(:,SLIMEpos(i)) > 0;
    metNames  = model.metNames(metPos);
    metStoich = model.S(metPos,SLIMEpos(i));
    
    %Find lipid species:
    pos_i = [];
    for j = 1:length(lipidNames)
        lipidPos = contains(model.metNames(metPos),lipidNames{j});
        if sum(lipidPos) > 0
            pos_i = j;
        end
    end
    
    %Find chains produced:
    for j = 1:length(chains)
        chainPos = strcmp(metNames,chains{j});
        if sum(chainPos) > 0
            composition(pos_i,j) = composition(pos_i,j) + flux*metStoich(chainPos)/mu;
        end
    end
end

%Find variability:
variability.min = zeros(size(composition));
variability.max = zeros(size(composition));
for i = 1:length(lipidNames)
    for j = 1:length(chains)
        [minVal,maxVal]      = lipidFVA(model_cons,lipidNames{i},chains{j});
        variability.min(i,j) = minVal/mu;
        variability.max(i,j) = maxVal/mu;
        disp(['Computing composition and variability: ' lipidNames{i} ' - ' chains{j}])
    end
end

%Generate output:
data.comp        = composition*1000;        %mg/gDW
data.var.min     = variability.min*1000;    %mg/gDW
data.var.max     = variability.max*1000;    %mg/gDW
data.vgluc       = sol.x(posGluc);          %1/h
data.mu          = mu;                      %1/h
data.maintenance = sol.x(posMaint);         %mmol/gDWh
data.netATP      = sol.x(posMaint)/mu;      %mmol/gDW

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%