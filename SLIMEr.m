%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = SLIMEr(model,abundances)
%
% Benjamín J. Sánchez. Last update: 2017-11-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = SLIMEr(model,data)

%Add new metabolites:
for i = 1:length(data.chainData.metNames)
    model = addLipidSpecies(model,data.chainData.metNames{i});
end
model = addLipidSpecies(model,'lipid - backbones');
model = addLipidSpecies(model,'lipid - tails');

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
model = changeChainComp(model,data.chainData);

%Change overall lipid pseudo-rxn:
tailID  = model.mets{strcmp(model.metNames,'lipid - tails [cytoplasm]')};
backID  = model.mets{strcmp(model.metNames,'lipid - backbones [cytoplasm]')};
lipidID = 's_1096[c]';
rxnName = 'lipid pseudoreaction - merge';
model   = addReaction(model, ...                    %model
                      {'r_2108',rxnName}, ...       %rxn
                      {tailID,backID,lipidID}, ...  %metabolites
                      [-1,    -1,    +1], ...       %stoichiometry
                      false, ...                    %reversibility
                      0, ...                        %LB
                      1000, ...                     %UB
                      0);                           %c

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%