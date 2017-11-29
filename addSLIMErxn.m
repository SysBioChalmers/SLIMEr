%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addSLIMErxn(model,met)
%
% Benjamín J. Sánchez. Last update: 2017-11-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = addSLIMErxn(model,rxnID)

rxnPos = strcmp(model.rxns,rxnID);

%Specific lipid:
specPos  = find(model.S(:,rxnPos) < 0);
specID   = model.mets{specPos};
specName = model.metNames{specPos};
rxnName  = [specName ' SLIME rxn'];
specName = specName(1:strfind(specName,'[')-2);

%Generic lipid (backbone):
backPos  = find(model.S(:,rxnPos) > 0);
backID   = model.mets{backPos};
backName = model.metNames{backPos};
backName = backName(1:strfind(backName,'[')-2);

%Define number of tails to add:
switch backName
    %Cases with format '(CXX)':
    case {'inositol-P-ceramide'
          'inositol phosphomannosylinositol phosphoceramide'
          'mannosylinositol phosphorylceramide'}
        tailsRxn = {};  %?????
        
    %Cases with specific erg-ester name:
    case 'ergosterol ester'
        esters = {'ergosteryl palmitoleate'     'C16:1'
                  'ergosteryl oleate'           'C18:1'};
        tailsRxn = esters(strcmp(esters(:,1),specName),2);
        
    %Cases with specific FA name:
    %OBS: For now, only marked ones are being added in the model
    case 'fatty acid'
        FAs = {'myristate'          'C14:0'     %
               'myristoleate'       'C14:1'
               'pentadecanoate'     'C15:0'
               'palmitate'          'C16:0'     %
               'palmitoleate'       'C16:1'     %
               'stearate'           'C18:0'     %
               'oleate'             'C18:1'     %
               'nonadecanoate'      'C19:0'
               'eicosanoate'        'C20:0'};
        tailsRxn = FAs(strcmp(FAs(:,1),specName),2);
        
    %Cases with format '(1-XX:Y, 2-XX:Y, ...)':
    case {'1-phosphatidyl-1D-myo-inositol'
          'phosphatidyl-L-serine'
          'phosphatidylcholine'
          'phosphatidylethanolamine'
          'triglyceride'}
        %Find all tails:
        tailPos = strfind(specName,':');
        tailsRxn   = cell(size(tailPos));
        for i = 1:length(tailPos)
            tailsRxn{i} = ['C' specName(tailPos(i)-2:tailPos(i)+1)];
        end
        
    %No tails added (for dolichol & complex sphingolipid):    
    otherwise
        tailsRxn = {};
end

%Find tail metabolites in model:
tailPos    = find(~cellfun(@isempty,strfind(model.metNames,' chain [cytoplasm]')));
tailIDs    = model.mets(tailPos)';
tailsModel = model.metNames(tailPos)';
tailCoeffs = zeros(size(tailIDs));

%Match to corresponding tail:
for i = 1:length(tailsRxn)
    tailName   = [tailsRxn{i} ' chain [cytoplasm]'];
    tailMatch  = strcmp(tailsModel,tailName);
    tailCoeffs = tailCoeffs + tailMatch;
end

%Create SLIME rxn (with same ID as previous ISA rxn but different name):
model = addReaction(model, ...                    %model
                    {rxnID,rxnName}, ...          %rxn
                    [specID,backID,tailIDs], ...  %metabolites
                    [-1,+1,tailCoeffs], ...       %stoichiometry
                    false, ...                    %reversibility
                    0, ...                        %LB
                    1000, ...                     %UB
                    0);                           %c

%Remove wrongly created field:
model = rmfield(model,'grRules');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%