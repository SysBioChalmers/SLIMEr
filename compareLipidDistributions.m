%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareLipidDistributions
%
% Benjamín J. Sánchez. Last update: 2017-11-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Original model:
model_original = load('yeast_7.8.mat');
model_original = model_original.model;

%Lipid data:
fid = fopen('lipid_data.csv');
lipidData = textscan(fid,'%s %s %s %f32','Delimiter',',','HeaderLines',1);
data.lipidData.metIDs    = lipidData{3};
data.lipidData.abundance = lipidData{4};
fclose(fid);

%Chain data:
fid = fopen('chain_data.csv');
chainData = textscan(fid,'%s %s %f32','Delimiter',',','HeaderLines',1);
data.chainData.metNames  = chainData{1};
data.chainData.formulas  = chainData{2};
data.chainData.abundance = chainData{3};
fclose(fid);

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