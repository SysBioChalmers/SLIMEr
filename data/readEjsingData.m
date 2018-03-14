%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = readEjsingData(i,model)
%
% Benjamín J. Sánchez. Last update: 2018-03-07
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

%Map each lipid to corresponding positions in the model. If at least one found, 
%get MW and transform units to abundance (assuming 8% of total lipid) if
%none found, filter out:
MWs = zeros(size(data.metNames));
for i = 1:length(data.metNames)
    pos = matchToModel(model,data.metNames{i});
    if sum(pos) > 0
        MWi    = getMWfromFormula(model.metFormulas(pos));
        MWs(i) = mean(MWi);
    end
end
data.metNames  = data.metNames(MWs > 0);
data.abundance = data.abundance(MWs > 0);
data.std       = data.std(MWs > 0);
MWs            = MWs(MWs > 0);
lipidContent   = 0.08;                                  %g(tot lipid)/gDW
lipidContent   = lipidContent*sum(data.abundance);      %g(tot lipid)/gDW, corrected
data.std       = data.std.*MWs;                         %g(lipid i)/mol(tot lipid)
data.abundance = data.abundance.*MWs;                  	%g(lipid i)/mol(tot lipid)
data.std       = data.std/sum(data.abundance);          %g(lipid i)/g(tot lipid)
data.abundance = data.abundance/sum(data.abundance);	%g(lipid i)/g(tot lipid)
data.std       = data.std*lipidContent;                 %g(lipid i)/gDW
data.abundance = data.abundance*lipidContent;           %g(lipid i)/gDW

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%