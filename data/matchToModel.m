%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pos = matchToModel(model,metName)
%
% Benjamín J. Sánchez. Last update: 2018-03-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pos = matchToModel(model,metName)

%Translating codes to names in model:
codes = {'TAG'          'triglyceride'                                      'endoplasmic reticulum membrane'
         'DAG'          'diglyceride'                                       'endoplasmic reticulum membrane'
         'PA'           'phosphatidate'                                     'endoplasmic reticulum membrane'
         'PS'           'phosphatidyl-L-serine'                             'endoplasmic reticulum membrane'
         'PE'           'phosphatidylethanolamine'                          'endoplasmic reticulum membrane'
         'PC'           'phosphatidylcholine'                               'endoplasmic reticulum membrane'
         'LPI'          'sn-2-acyl-1-lysophosphatidylinositol'              'endoplasmic reticulum membrane'
         'PI'           '1-phosphatidyl-1D-myo-inositol'                    'cytoplasm'
         'PG'           'phosphatidylglycerol'                              'mitochondrial membrane'
         'CL'           'cardiolipin'                                       'mitochondrial membrane'
         'LCB 18:0;2'	'sphinganine'                                       'endoplasmic reticulum'
         'LCB 18:0;3'   'phytosphingosine'                                  'endoplasmic reticulum'
         'LCBP 18:0;2'	'sphinganine 1-phosphate'                           'endoplasmic reticulum'
         'LCBP 18:0;3'  'phytosphingosine 1-phosphate'                      'endoplasmic reticulum'
         'Cer'          'ceramide'                                          'Golgi'
         'IPC'          'inositol-P-ceramide'                               'Golgi'
         'MIPC'         'mannosylinositol phosphorylceramide'               'Golgi'
         'M(IP)2C'      'inositol phosphomannosylinositol phosphoceramide'  'Golgi'
         'Ergosterol'   'ergosterol'                                        'cytoplasm'};

pos = false(size(model.mets));
if contains(metName,'LCB') || contains(metName,'Ergost')
    %Direct match (LCB 18:0;2, LCB 18:0;3, LCBP 18:0;2, LCBP 18:0;3, Ergosterol):
    codesPos = strcmp(codes(:,1),metName);
    metName  = [codes{codesPos,2} ' [' codes{codesPos,3} ']'];
    pos      = strcmp(model.metNames,metName);
else
    %Split in backbone and tails:
    parts    = strsplit(metName,' ');
    backCode = parts{1};
    tailCode = parts{2};
        
    %Get full backbone and compartment name:
    codesPos = strcmp(codes(:,1),backCode);
    backName = codes{codesPos,2};
    compName = codes{codesPos,3};
        
    %Find tails and act accordingly:
    tails = strsplit(tailCode,'-');
    switch backCode
        case {'TAG', 'DAG', 'PA', 'PS', 'PE', 'PC', 'LPI', 'PI', 'PG', 'CL'}
            %Update pos with all possible combinations:
            tailOptions = perms(tails);
            [m,n]       = size(tailOptions);
            for i = 1:m
                tail_i      = strcat(num2str((1:n)'),'-',tailOptions(i,:)');
                metName     = [backName ' (' strjoin(tail_i,', ') ') [' compName ']'];
                metPos      = strcmp(model.metNames,metName);
                pos(metPos) = true;
            end
        case {'Cer', 'IPC', 'MIPC', 'M(IP)2C'}
            %Match to the corresponding variant:
            sphingos = {'1'    'A'     '2'     '0'     %cer1  = cerA  = ;2-;0
                        '2'    'B'     '3'     '0'     %cer2  = cerB  = ;3-;0
                        '3'    'C'     '3'     '1'     %cer3  = cerC  = ;3-;1
                        '4'    'D'     '3'     '2'     %cer4  = cerD  = ;3-;2
                        '2'''  'B'''   '2'     '1'};   %cer2' = cerB' = ;2-;1
            shortTail = tails{1}(end);
            longTail  = tails{2}(end);
            posSphing = strcmp(sphingos(:,3),shortTail).*strcmp(sphingos(:,4),longTail) == 1;
            metName   = [backName '-' sphingos{posSphing,1} ' (C' tails{2}(1:2) ') [' compName ']'];    %Cer
            if sum(strcmp(model.metNames,metName)) == 0
                metName = [backName ' ' sphingos{posSphing,2} ' (C' tails{2}(1:2) ') [' compName ']'];  %IPC, MIPC & M(IP)2C
            end
            pos = strcmp(model.metNames,metName);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%