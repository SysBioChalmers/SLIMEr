%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = convertEjsingData(data,model)
%
% Benjamín J. Sánchez. Last update: 2018-03-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = convertEjsingData(data,model)

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
data.abundance = data.abundance(MWs > 0)/100;           %mol(lipid i)/mol(tot lipid)
data.std       = data.std(MWs > 0)/100;                 %mol(lipid i)/mol(tot lipid)
MWs            = MWs(MWs > 0);
lipidContent   = 0.08;                                  %g(tot lipid)/gDW
lipidContent   = lipidContent*sum(data.abundance);      %g(tot lipid)/gDW, corrected
data.std       = data.std.*MWs;                         %g(lipid i)/mol(tot lipid)
data.abundance = data.abundance.*MWs;                  	%g(lipid i)/mol(tot lipid)
data.std       = data.std/sum(data.abundance);          %g(lipid i)/g(tot lipid)
data.abundance = data.abundance/sum(data.abundance);	%g(lipid i)/g(tot lipid)
data.std       = data.std*lipidContent;                 %g(lipid i)/gDW, corrected
data.abundance = data.abundance*lipidContent;           %g(lipid i)/gDW, corrected

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%