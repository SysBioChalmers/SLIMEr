%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createNewModels
%
% Benjamín J. Sánchez. Last update: 2018-01-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
addpath('../simulations')
addpath('../data')

%Original model:
model_original = load('yeast_7.8.mat');
model_original = model_original.model;

%Read data:
data = readData;

%Model with lipid composition corrected:
model_correctedComp = SLIMEr(model_original,data,false);

%Model with both lipid and chain length constrained to data:
model_SLIMEr = SLIMEr(model_original,data,true);

%Make abundances be consistent:
[model_SLIMEr,k]    = scaleAbundancesInModel(model_SLIMEr,data);
model_correctedComp = adjustModel(model_correctedComp,k,false);

%Correct the rest of the composition to be consistent:
model_correctedComp = changeOtherComp(model_correctedComp,data.otherData);
model_SLIMEr        = changeOtherComp(model_SLIMEr,data.otherData);
rmpath('../simulations')
rmpath('../data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%