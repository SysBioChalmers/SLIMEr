%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = SLIMEr(model,data,includeTails)
%
% Benjamín J. Sánchez. Last update: 2018-03-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = SLIMEr(model,data,includeTails)

%Add new metabolites:
for i = 1:length(data.chainData.metNames)
    metName = [data.chainData.metNames{i} ' [cytoplasm]'];
    model   = addLipidSpecies(model,metName,data.chainData.formulas{i},~includeTails);
end

%Create lipid pseudo-rxn for backbones:
model = addLipidSpecies(model,'lipid - backbones [cytoplasm]','NA',includeTails);
model = changeLipidComp(model,data.lipidData);

%Create lipid pseudo-rxn for tails:
if includeTails
    model = addLipidSpecies(model,'lipid - tails [cytoplasm]','NA',includeTails);
    model = changeChainComp(model,data.chainData);
end

%Add SLIME reactions replacing existing ISA rxns:
rxnIDs   = model.rxns;
rxnNames = model.rxnNames;
toDelete = false(size(rxnIDs));
for i = 1:length(rxnIDs)
    if contains(rxnNames{i},'isa ')
        [model,toDelete(i)] = addSLIMErxn(model,rxnIDs{i},[]);
    end
end

%Delete ISA rxns not used:
model = removeRxns(model,rxnIDs(toDelete));

%Add new SLIME reactions (with no corresponding ISA rxns):
metIDs   = model.mets;
metNames = model.metNames;
for i = 1:length(metIDs)
    backName = getBackboneName(metNames{i});
    if ~isempty(backName)
        [model,~] = addSLIMErxn(model,[],metIDs{i});
    end
end

%Change overall lipid pseudo-rxn:
backID  = model.mets{strcmp(model.metNames,'lipid - backbones [cytoplasm]')};
lipidID = 's_1096[c]';
if includeTails
    tailID = model.mets{strcmp(model.metNames,'lipid - tails [cytoplasm]')};
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

%Remove GAM requirement:
GAM    = 0;
bioRxn = strcmp(model.rxnNames,'yeast 8 biomass pseudoreaction');
ATPpos = strcmp(model.metNames,'ATP [cytoplasm]');
H2Opos = strcmp(model.metNames,'H2O [cytoplasm]');
ADPpos = strcmp(model.metNames,'ADP [cytoplasm]');
Hpos   = strcmp(model.metNames,'H+ [cytoplasm]');
Ppos   = strcmp(model.metNames,'phosphate [cytoplasm]');
model.S(ATPpos,bioRxn) = -GAM;
model.S(H2Opos,bioRxn) = -GAM;
model.S(ADPpos,bioRxn) = +GAM;
model.S(Hpos,bioRxn)   = +GAM;
model.S(Ppos,bioRxn)   = +GAM;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%