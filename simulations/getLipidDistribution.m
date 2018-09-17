%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = getLipidDistribution(model,lipidNames,chains,fluxData,getFullVar)
%
% Benjamin J. Sanchez. Last update: 2018-09-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = getLipidDistribution(model,lipidNames,chains,fluxData,getFullVar)

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

%Compute full FVA:
if getFullVar
    fullVar = zeros(size(model_cons.rxns));
    for i = 1:length(fullVar)
        %Change objective function and perform FVA:
        model_i = changeObjective(model_cons,model_cons.rxns{i},+1);
        sol_min = optimizeCbModel(model_i,'min');
        sol_max = optimizeCbModel(model_i,'max');
        fullVar(i) = sol_max.x(i) - sol_min.x(i);
        if rem(i,50) == 0
            disp(['Computing FVA for reference: ' num2str(i) '/' ...
                num2str(length(model.rxns)) ' rxns done'])
        end
    end
    data.fullVar = fullVar;  %mmol/gDWh
end

%Find lipid variability:
lipVar.min = zeros(size(composition));
lipVar.max = zeros(size(composition));
for i = 1:length(lipidNames)
    for j = 1:length(chains)
        [minVal,maxVal] = lipidFVA(model_cons,lipidNames{i},chains{j});
        lipVar.min(i,j) = minVal/mu;
        lipVar.max(i,j) = maxVal/mu;
    end
end

%Generate output:
data.comp        = composition*1000;	%mg/gDW
data.lipVar.min  = lipVar.min*1000;     %mg/gDW
data.lipVar.max  = lipVar.max*1000;     %mg/gDW
data.vgluc       = sol.x(posGluc);      %1/h
data.mu          = mu;                  %1/h
data.maintenance = sol.x(posMaint);     %mmol/gDWh
data.netATP      = sol.x(posMaint)/mu;	%mmol/gDW

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
