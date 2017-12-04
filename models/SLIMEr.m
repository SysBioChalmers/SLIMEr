%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = SLIMEr(model,data,includeTails)
%
% Benjamín J. Sánchez. Last update: 2017-12-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = SLIMEr(model,data,includeTails)

%Add new metabolites:
for i = 1:length(data.chainData.metNames)
    model = addLipidSpecies(model,data.chainData.metNames{i},data.chainData.formulas{i},~includeTails);
end
model = addLipidSpecies(model,'lipid - backbones','NA',includeTails);
model = addLipidSpecies(model,'lipid - tails','NA',includeTails);

%Add SLIME reactions:
for i = 1:length(model.rxns)
    if ~isempty(strfind(model.rxnNames{i},'isa '))
        model = addSLIMErxn(model,model.rxns{i});
        printRxnFormula(model,model.rxns{i},true,true,true);
    end
end

%Create lipid pseudo-rxn for backbones:
model = changeLipidComp(model,data.lipidData);

%Create lipid pseudo-rxn for tails:
if includeTails
    model  = changeChainComp(model,data.chainData);
end

%Change overall lipid pseudo-rxn:
tailID  = model.mets{strcmp(model.metNames,'lipid - tails [cytoplasm]')};
backID  = model.mets{strcmp(model.metNames,'lipid - backbones [cytoplasm]')};
lipidID = 's_1096[c]';
if includeTails
    mets   = {tailID,backID,lipidID};
    coeffs = [-1,    -1,    +1];
else
    mets   = {backID,lipidID};
    coeffs = [-1,    +1];
end
rxnName = 'lipid pseudoreaction - merge';
model   = addReaction(model, ...                    %model
                      {'r_2108',rxnName}, ...       %rxn
                      mets, ...                     %metabolites
                      coeffs, ...                   %stoichiometry
                      false, ...                    %reversibility
                      0, ...                        %LB
                      1000, ...                     %UB
                      0);                           %c

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%