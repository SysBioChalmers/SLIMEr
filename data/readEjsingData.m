%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = readEjsingData(i,model)
%
% Benjamín J. Sánchez. Last update: 2018-03-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = readEjsingData(i,model)

%Load data:
fid            = fopen('fullData_Ejsing2009.csv');
lipidData      = textscan(fid,['%s ' repmat('%f32 ',[1,15]) '%f32'],'Delimiter',',','HeaderLines',1);
data.metNames  = lipidData{1};
data.abundance = lipidData{2*i};
data.std       = lipidData{2*i+1};
fclose(fid);

%Filter out:
%  * any species with 10, 12, 14, 20 or 22 carbons
%  * any double saturation (":2")
%  * lipid classes missing in model (LPA, LPS, LPE, LPC, LIPC, ergostadienol)
filterOut = contains(data.metNames,'10')   + contains(data.metNames,'12') + ...
            contains(data.metNames,'14')   + contains(data.metNames,'20') + ...
            contains(data.metNames,'22')   + contains(data.metNames,':2') + ...
            contains(data.metNames,'LPA')  + contains(data.metNames,'LPS') + ...
            contains(data.metNames,'LPE')  + contains(data.metNames,'LPC') + ...
            contains(data.metNames,'LIPC') + contains(data.metNames,'Ergostadienol');
data.metNames  = data.metNames(~filterOut);
data.abundance = data.abundance(~filterOut);
data.std       = data.std(~filterOut);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%