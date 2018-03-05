%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model,toDelete] = addSLIMErxn(model,rxnID,specID)
%
% Benjamín J. Sánchez. Last update: 2018-03-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,toDelete] = addSLIMErxn(model,rxnID,specID)

%Already existing ISA-rxns:
if isempty(specID)
    rxnPos  = strcmp(model.rxns,rxnID);
    printRxnFormula(model,model.rxns{rxnPos},true,true,true);
    specPos = find(model.S(:,rxnPos) < 0);
    specID  = model.mets{specPos};
end

%Non existing ISA rxn:
if isempty(rxnID)
    rxnID   = ['r_' getNewIndex(model.rxns)];
    specPos = strcmp(model.mets,specID);
end

%Specific lipid:
specName = model.metNames{specPos};
rxnName  = [specName ' SLIME rxn'];

%Generic lipid (backbone):
backName = getBackboneName(specName);
backPos  = strcmp(model.metNames,backName);
toDelete = false;
if isempty(backName)
    disp('removing')    %dolichol and any other weird ISA rxn
    toDelete = true;
    return
elseif sum(backPos) > 0
    backID = model.mets{backPos};       %backbone already exists
else
    model   = addLipidSpecies(model,backName,'',false); %add formula later
    backPos = strcmp(model.metNames,backName);
    backID  = model.mets{backPos};
end
backName = backName(1:strfind(backName,'[')-2);
specName = specName(1:strfind(specName,'[')-2);

%Stoich. coeff. for backbone: molecular weight of specific species
specMW = getMWfromFormula(model.metFormulas(specPos));

%Define number of tails to add:
switch backName       
    %Cases with format '(1-XX:Y, 2-XX:Y, ...)':
    case {'1-phosphatidyl-1D-myo-inositol'                      %PI
          'sn-2-acyl-1-lysophosphatidylinositol'                %LPI
          'phosphatidyl-L-serine'                               %PS
          'phosphatidylcholine'                                 %PC
          'phosphatidylethanolamine'                            %PE
          'phosphatidate'                                       %PA
          'phosphatidylglycerol'                                %PG
          'cardiolipin'                                         %CL
          'diglyceride'                                         %DAG
          'triglyceride'}                                       %TAG
        %Find all tails:
        tailPos  = strfind(specName,':');
        tailsRxn = cell(size(tailPos));
        for i = 1:length(tailPos)
            tailsRxn{i} = ['C' specName(tailPos(i)-2:tailPos(i)+1)];
        end
        
    %Cases with format '(CXX)', for very long fatty acid chains (24 & 26 carbons):
    case {'ceramide'                                            %Cer
          'inositol-P-ceramide'                                 %IPC
          'inositol phosphomannosylinositol phosphoceramide'    %MIPC
          'mannosylinositol phosphorylceramide'}                %M(IP)2C
        if contains(specName,'(C24)')
            tailsRxn = {'C18:0';'C24:0'};
        elseif contains(specName,'(C26)')
            tailsRxn = {'C18:0';'C26:0'};
        end
        
    %Cases with specific names:
    case {'fatty acid'
          'ergosterol ester'
          'long-chain base'
          'long-chain base phosphate'}
        species = {'myristate'                      'C14:0'     %FA
                   'myristoleate'                   'C14:1'     %FA
                   'pentadecanoate'                 'C15:0'     %FA
                   'palmitate'                      'C16:0'     %FA
                   'palmitoleate'                   'C16:1'     %FA
                   'stearate'                       'C18:0'     %FA
                   'oleate'                         'C18:1'     %FA
                   'nonadecanoate'                  'C19:0'     %FA
                   'eicosanoate'                    'C20:0'     %FA
                   'ergosteryl palmitoleate'        'C16:1'     %erg-ester
                   'ergosteryl oleate'              'C18:1'     %erg-ester
                   'sphinganine'                    'C18:0'     %LCB
                   'phytosphingosine'               'C18:0'     %LCB
                   'sphinganine 1-phosphate'        'C18:0'     %LCBP
                   'phytosphingosine 1-phosphate'	'C18:0'};   %LCBP
        tailsRxn = species(strcmp(species(:,1),specName),2);
        
    %No tails added (for dolichol & "complex sphingolipid"):    
    otherwise
        error('smth went wrong')
end

%Find tail metabolites in model:
tailPos    = find(~cellfun(@isempty,strfind(model.metNames,' chain [cytoplasm]')));
tailIDs    = model.mets(tailPos)';
tailsModel = model.metNames(tailPos)';
tailsMWs   = getMWfromFormula(model.metFormulas(tailPos)');
tailCoeffs = zeros(size(tailIDs));

%Match to corresponding tail:
for i = 1:length(tailsRxn)
    tailName   = [tailsRxn{i} ' chain [cytoplasm]'];
    tailMatch  = strcmp(tailsModel,tailName);
    tailCoeffs = tailCoeffs + tailMatch.*tailsMWs;
end

%Create SLIME rxn (with same ID as previous ISA rxn but different name):
model = addReaction(model,rxnID, ...
                    'reactionName', rxnName, ...
                    'metaboliteList', [specID,backID,tailIDs], ...
                    'stoichCoeffList', [-1,specMW,tailCoeffs], ...
                    'reversible', false, ...
                    'lowerBound', 0, ...
                    'upperBound', 1000, ...
                    'checkDuplicate', true);

%Remove wrongly created field:
try
    model = rmfield(model,'grRules');
    printRxnFormula(model,rxnID,true,true,true);
catch
    disp(['Repeated: ' rxnName])
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%