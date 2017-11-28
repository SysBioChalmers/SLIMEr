%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareLipidDistributions
%
% Benjamín J. Sánchez. Last update: 2017-11-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Original model:
model_original = load('yeast_7.7.mat');
model_original = model_original.model;

%Data:
rxnPos = strcmp(model_original.rxns,'r_2108');
data.lipidData.metIDs    = model_original.mets(model_original.S(:,rxnPos) < 0);            %format: s_XXXX[abc]
data.lipidData.abundance = abs(model_original.S(model_original.S(:,rxnPos) < 0,rxnPos));   %mmol/gDW (how??)
data.chainData.metNames  = {'C16:0 chain'
                            'C16:1 chain'
                            'C18:0 chain'
                            'C18:1 chain'};     %format: CXX:Y chain
data.chainData.abundance = [0.005/256;0.005/254;0.005/284;0.005/282]*1000;   %5 mg/gDW for each species = 2% in mass for all lipids

%Model with lipid composition corrected:
model_correctedComp = changeLipidComp(model_original,data.lipidData);

%Model with both lipid and chain length constrained to data:
model_SLIMEr = SLIMEr(model_original,data);

%Simulate models:
sol_original      = optimizeCbModel(model_original);
sol_correctedComp = optimizeCbModel(model_correctedComp);
sol_SLIMEr        = optimizeCbModel(model_SLIMEr);

%Compare distributions:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%