%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createNewModels
%
% Benjamín J. Sánchez. Last update: 2018-01-18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
addpath('../simulations')
addpath('../data')

%Original model:
model_original = load('yeast_7.8.mat');
model_original = model_original.model;

%Create 2 models for each of the 10 conditions:
model_correctedComp = cell(1,10);
model_SLIMEr        = cell(1,10);
GAMpol              = zeros(1,10);
k                   = zeros(1,10);
for i = 1:10
    %Read data:
    data = readData(i);
    
    %Model with lipid composition corrected:
    model_correctedComp{i} = SLIMEr(model_original,data,false);
    
    %Model with both lipid and chain length constrained to data:
    model_SLIMEr{i} = SLIMEr(model_original,data,true);
    
    %Make abundances be consistent:
    [model_SLIMEr{i},k(i)] = scaleAbundancesInModel(model_SLIMEr{i},data);
    model_correctedComp{i} = adjustModel(model_correctedComp{i},k(i),false);
    
    %Correct the rest of the composition to be consistent:
    [model_correctedComp{i},~]  = changeOtherComp(model_correctedComp{i},data);
    [model_SLIMEr{i},GAMpol(i)] = changeOtherComp(model_SLIMEr{i},data);
end

rmpath('../simulations')
rmpath('../data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%