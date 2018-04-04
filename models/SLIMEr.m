%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = SLIMEr(model,data,includeTails)
%
% Benjamín J. Sánchez. Last update: 2018-03-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = SLIMEr(model,data,includeTails)

%Add new backbones:
metIDs   = model.mets;
metNames = model.metNames;
for i = 1:length(metIDs)
    backName = getBackboneName(metNames{i});
    if ~isempty(backName)
        if ~ismember(backName,model.metNames)
            model = addLipidSpecies(model,backName,'',false);
        end
        %Add transport rxn to cytoplasm for non cytoplasmic backbones:
        if ~contains(backName,'cytoplasm')
            cytoName = [backName(1:strfind(backName,'[')) 'cytoplasm]'];
            if ~ismember(cytoName,model.metNames)
                model = addLipidSpecies(model,cytoName,'',false);
            end
            backID    = model.mets(strcmp(model.metNames,backName));
            cytoID    = model.mets(strcmp(model.metNames,cytoName));
            transID   = ['r_' getNewIndex(model.rxns)];
            transName = [backName(1:strfind(backName,'[')-1) 'transport'];
            model     = addReaction(model,transID, ...
                                    'reactionName', transName, ...
                                    'metaboliteList', [backID,cytoID], ...
                                    'stoichCoeffList', [-1,+1], ...
                                    'reversible', false, ...
                                    'lowerBound', 0, ...
                                    'upperBound', 1000, ...
                                    'checkDuplicate', true);
            printRxnFormula(model,transID,true,true,true);
        end
    end
end

%Get backbone IDs:
metIDs = cell(size(data.lipidData.metNames));
for i = 1:length(data.lipidData.metNames)
    metName   = [data.lipidData.metNames{i} ' [cytoplasm]'];
    metIDs{i} = model.mets{strcmp(model.metNames,metName)};
end
data.lipidData.metIDs = metIDs;

%Add chains:
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
    elseif strcmp(rxnNames{i},'complex sphingolipid transport')
        toDelete(i) = true;
    end
end

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
model   = addReaction(model, 'r_2108', ...
                      'reactionName', rxnName, ...
                      'metaboliteList', mets, ...
                      'stoichCoeffList', coeffs, ...
                      'reversible', false, ...
                      'lowerBound', 0, ...
                      'upperBound', 1000);
                  
%Delete ISA rxns not used:
model = removeRxns(model,rxnIDs(toDelete));

%Remove GAM requirement:
GAM    = 0;
bioRxn = strcmp(model.rxnNames,'biomass pseudoreaction');
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

%Remove wrongly created field:
model = rmfield(model,'grRules');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%