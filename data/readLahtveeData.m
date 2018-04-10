%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = readLahtveeData(i)
%
% Benjamín J. Sánchez. Last update: 2018-04-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = readLahtveeData(i)

%Lipid data:
fid = fopen('lipidData_Lahtvee2016.csv');
lipidData = textscan(fid,[repmat('%s ',[1,3]) repmat('%f32 ',[1,9]) '%f32'],'Delimiter',',','HeaderLines',1);
data.lipidData.metAbbrev = lipidData{1};
data.lipidData.metNames  = lipidData{2};
data.lipidData.metIDs    = lipidData{3};
data.lipidData.abundance = lipidData{3+i};
fclose(fid);

%Chain data:
fid = fopen('chainData_Lahtvee2016.csv');
chainData = textscan(fid,[repmat('%s ',[1,2]) repmat('%f32 ',[1,19]) '%f32'],'Delimiter',',','HeaderLines',1);
data.chainData.metNames  = chainData{1};
data.chainData.formulas  = chainData{2};
data.chainData.abundance = chainData{1+2*i};
data.chainData.std       = chainData{2+2*i};
fclose(fid);

%Other composition data:
fid = fopen('compData_Lahtvee2016.csv');
otherData = textscan(fid,[repmat('%s ',[1,2]) repmat('%f32 ',[1,9]) '%f32'],'Delimiter',',','HeaderLines',1);
data.otherData.metIDs    = otherData{2};
data.otherData.abundance = otherData{2+i};
fclose(fid);

%Flux data:
fid = fopen('fluxData_Lahtvee2016.csv');
fluxData = textscan(fid,[repmat('%s ',[1,2]) repmat('%f32 ',[1,19]) '%f32'],'Delimiter',',','HeaderLines',1);
data.fluxData.rxnIDs   = fluxData{2};
data.fluxData.averages = fluxData{2*i+1};
data.fluxData.stdevs   = fluxData{2*i+2};
fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%