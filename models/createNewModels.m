%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createNewModels
%
% Benjamín J. Sánchez. Last update: 2018-01-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
addpath('../simulations')

%Original model:
model_original = load('yeast_7.8.mat');
model_original = model_original.model;

%Lipid data:
fid = fopen('../data/lipidData_Lahtvee2016.csv');
lipidData = textscan(fid,'%s %s %s %f32','Delimiter',',','HeaderLines',1);
data.lipidData.metIDs    = lipidData{3};
data.lipidData.abundance = lipidData{4};
fclose(fid);

%Chain data:
fid = fopen('../data/chainData_Lahtvee2016.csv');
chainData = textscan(fid,'%s %s %f32','Delimiter',',','HeaderLines',1);
data.chainData.metNames  = chainData{1};
data.chainData.formulas  = chainData{2};
data.chainData.abundance = chainData{3};
fclose(fid);

%Other composition data:
fid = fopen('../data/compData_Lahtvee2016.csv');
otherData = textscan(fid,'%s %s %f32','Delimiter',',','HeaderLines',1);
data.otherData.metIDs    = otherData{2};
data.otherData.abundance = otherData{3};
fclose(fid);

%Flux data:
fid = fopen('../data/fluxData_Lahtvee2016.csv');
fluxData = textscan(fid,'%s %s %f32 %f32','Delimiter',',','HeaderLines',1);
data.fluxData.rxnIDs   = fluxData{2};
data.fluxData.averages = fluxData{3};
data.fluxData.stdevs   = fluxData{4};
fclose(fid);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%