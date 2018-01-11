%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% readData
%
% Benjamín J. Sánchez. Last update: 2018-01-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = readData

%Lipid data:
fid = fopen('lipidData_Lahtvee2016.csv');
lipidData = textscan(fid,'%s %s %s %f32','Delimiter',',','HeaderLines',1);
data.lipidData.metAbbrev = lipidData{1};
data.lipidData.metNames  = lipidData{2};
data.lipidData.metIDs    = lipidData{3};
data.lipidData.abundance = lipidData{4};
fclose(fid);

%Chain data:
fid = fopen('chainData_Lahtvee2016.csv');
chainData = textscan(fid,'%s %s %f32','Delimiter',',','HeaderLines',1);
data.chainData.metNames  = chainData{1};
data.chainData.formulas  = chainData{2};
data.chainData.abundance = chainData{3};
fclose(fid);

%Other composition data:
fid = fopen('compData_Lahtvee2016.csv');
otherData = textscan(fid,'%s %s %f32','Delimiter',',','HeaderLines',1);
data.otherData.metIDs    = otherData{2};
data.otherData.abundance = otherData{3};
fclose(fid);

%Flux data:
fid = fopen('fluxData_Lahtvee2016.csv');
fluxData = textscan(fid,'%s %s %f32 %f32','Delimiter',',','HeaderLines',1);
data.fluxData.rxnIDs   = fluxData{2};
data.fluxData.averages = fluxData{3};
data.fluxData.stdevs   = fluxData{4};
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%