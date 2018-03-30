%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createNewModels
%
% Benjamín J. Sánchez. Last update: 2018-03-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
addpath('../simulations')
addpath('../data')

%Start COBRA:
initCobraToolbox

%Original model:
model = load('yeastGEM.mat');
model = model.model;

%Create 2 models for each of the 10 conditions in the stress dataset:
model_corrComp = cell(1,10);
model_SLIMEr   = cell(1,10);
GAMpol         = zeros(1,10);
k              = zeros(1,10);
for i = 1:10
    data = readLahtveeData(i);
    [model_corrComp{i},model_SLIMEr{i},k(i),GAMpol(i)] = modelsFromData(model,data,'backbones');
end

%Create 2 model for each of the validation studies:
model_corrComp_val = cell(1,8);
model_SLIMEr_val   = cell(1,8);
for i = 1:8
    data = readEjsingData(i);
    data = convertEjsingData(data,model,true);
    [model_corrComp_val{i},model_SLIMEr_val{i},~,~] = modelsFromData(model,data,'tails');
end

rmpath('../simulations')
rmpath('../data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%