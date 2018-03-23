%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createNewModels
%
% Benjamín J. Sánchez. Last update: 2018-03-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
addpath('../simulations')
addpath('../data')

%Original model:
model = load('yeast_7.8.mat');
model = model.model;

%Create 2 models for each of the 10 conditions:
model_corrComp = cell(1,10);
model_SLIMEr        = cell(1,10);
GAMpol              = zeros(1,10);
k                   = zeros(1,10);
for i = 1:10
    data = readLahtveeData(i);
    [model_corrComp{i},model_SLIMEr{i},k(i),GAMpol(i)] = modelsFromData(model,data);
end

rmpath('../simulations')
rmpath('../data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%