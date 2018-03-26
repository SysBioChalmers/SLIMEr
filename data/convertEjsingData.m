%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = convertEjsingData(data,model,condense)
%
% Benjamín J. Sánchez. Last update: 2018-03-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = convertEjsingData(data_old,model,condense)

%Map each lipid to corresponding positions in the model:
MWs        = zeros(size(data_old.metNames));
backbones  = cell(size(data_old.metNames));
for i = 1:length(data_old.metNames)
    pos = matchToModel(model,data_old.metNames{i});
    if sum(pos) > 0
        MWi          = getMWfromFormula(model.metFormulas(pos));
        MWs(i)       = mean(MWi);
        backbone     = model.metNames(pos);
        backbones{i} = getBackboneName(backbone{1});
        backbones{i} = backbones{i}(1:strfind(backbones{i},'[')-2);
    end
end

%Filter out mets with no MW computed (i.e. not in model):
backbones          = backbones(MWs > 0);
data_old.metNames  = data_old.metNames(MWs > 0);
data_old.abundance = data_old.abundance(MWs > 0)/100;	%mol(lipid i)/mol(tot lipid)
data_old.std       = data_old.std(MWs > 0)/100;         %mol(lipid i)/mol(tot lipid)
MWs                = MWs(MWs > 0);

if condense
    %Backbones: First convert to g/gDW and then add up
    data.lipidData = changeUnits(data_old,MWs);
    data.lipidData = squashLipids(data.lipidData,backbones);
    
    %Tails: First get tail names & MWs
    fid       = fopen('chainData_Lahtvee2016.csv');
    chainData = textscan(fid,[repmat('%s ',[1,2]) repmat('%f32 ',[1,9]) '%f32'],'Delimiter',',','HeaderLines',1);
    fclose(fid);
    chainNames    = chainData{1};
    tempNames     = strrep(chainNames,'C','');
    tempNames     = strrep(tempNames,' chain','');
    chainFormulas = chainData{2};
    chainMWs      = getMWfromFormula(chainFormulas);
    %Now add up and correct fields:
    data.chainData          = squashLipids(data_old,tempNames);
    data.chainData.formulas = chainFormulas;
    data.chainData.metNames = chainNames;
    %Finally, convert to g/gDW:
    data.chainData = changeUnits(data.chainData,chainMWs);
else
    %Only change units to g/gDW:
    data = changeUnits(data_old,MWs);
end

%Add missing data: composition and fluxes
data.otherData.metIDs    = {'protein';'RNA'};
data.otherData.abundance = [0.5;0.06];              %Average values from literature
data.fluxData.rxnIDs     = {'r_1714';'r_2111'};     %Glucose & biomass
data.fluxData.averages   = [-20.4;0.41];            %Estimated by ecYeast7
data.fluxData.stdevs     = [20.4;0.41]*0.01;        %Assume 1% error in measurement

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = changeUnits(data,MWs)

%Transform units to abundance (assuming 8% of total lipid)
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

function data = squashLipids(data_old,metNames)

%Initialize variables:
data.metNames  = unique(metNames);
data.abundance = zeros(size(data.metNames));
data.std       = zeros(size(data.metNames));

%Determine if backbones or tails are being squashed:
if length(data.metNames) == length(metNames)
    isTail = true;
else
    isTail = false;
end

%Squash each species by adding all abundances and averaging std:
for i = 1:length(data.metNames)
    if isTail
        hits = zeros(size(data_old.metNames));
        for j = 1:length(data_old.metNames)
            hits(j) = length(strfind(data_old.metNames{j},data.metNames{i}));
        end
    else
        hits = strcmp(metNames,data.metNames{i});
    end
    data.abundance(i) = sum(data_old.abundance.*hits);
    data.std(i)       = mean(data_old.std.*hits);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%