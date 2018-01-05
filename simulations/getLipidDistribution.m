%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = getLipidDistribution(model,lipidNames,chains)
%
% Benjamín J. Sánchez. Last update: 2018-01-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = getLipidDistribution(model,lipidNames,chains)

%Simulate model:
sol = simulateGrowth(model,'growth');

%Find growth:
posG = strcmp(model.rxnNames,'D-glucose exchange');
posX = strcmp(model.rxnNames,'growth');
mu   = sol.x(posX);

%Find energy requirements:
[ATP,NADH,NADPH] = getEnergyexpenditure(model,sol.x/mu);

for i = 1:length(chains)
    chains{i} = ['C' chains{i} ' chain [cytoplasm]'];
end

composition = zeros(length(lipidNames),length(chains));

%Go through all SLIME rxns to find abundances:
SLIMEpos = find(~cellfun(@isempty,strfind(model.rxnNames,'SLIME rxn')));
for i = 1:length(SLIMEpos)
    %Find flux and all metabolites produced in each SLIME rxn:
    flux      = sol.x(SLIMEpos(i));
    metPos    = model.S(:,SLIMEpos(i)) > 0;
    metNames  = model.metNames(metPos);
    metStoich = model.S(metPos,SLIMEpos(i));
    
    %Find lipid species:
    pos_i = [];
    for j = 1:length(lipidNames)
        lipidPos = ~cellfun(@isempty,strfind(model.metNames(metPos),lipidNames{j}));
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

%Fix glucose and biomass:
model = changeRxnBounds(model,model.rxns(posG),-1,'b');
model = changeRxnBounds(model,model.rxns(posX),mu*0.999,'l');

%Find variability:
variability.min = zeros(size(composition));
variability.max = zeros(size(composition));
for i = 1:length(lipidNames)
    for j = 1:length(chains)
        [minVal,maxVal]      = lipidFVA(model,lipidNames{i},chains{j});
        variability.min(i,j) = minVal/mu;
        variability.max(i,j) = maxVal/mu;
        disp(['Computing composition and variability: ' lipidNames{i} ' - ' chains{j}])
    end
end

data.comp    = composition*1000;        %mg/gDW
data.var.min = variability.min*1000;    %mg/gDW
data.var.max = variability.max*1000;    %mg/gDW
data.mu      = mu;
data.ATP     = ATP;
data.NADH    = NADH;
data.NADPH   = NADPH;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%