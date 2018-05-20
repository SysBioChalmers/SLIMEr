%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stoich = getStoichFromFormula(metFormulas,element)
% Returns the stoichiometry of a given element in a set of chemical formulas
%
% Benjamín J. Sánchez. Last update: 2018-05-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stoich = getStoichFromFormula(metFormulas,element)

%Calculate the MW for each formula:
stoich = zeros(size(metFormulas));
for i = 1:length(stoich)
    formula = metFormulas{i};
    elePos  = strfind(formula,element);
    if length(elePos) > 1
        error(['Non-standard formula: ' formula])
    elseif isempty(elePos)
        stoich(i) = 0;
    elseif elePos == length(formula)
        stoich(i) = 1;
    else
        rest    = formula(elePos+1:end);
        nextPos = regexp(rest,'[A-Z]');
        if isempty(nextPos)
            stoich(i) = str2double(rest);
        elseif nextPos(1) == 1
            stoich(i) = 1;
        else
            stoich(i) = str2double(rest(1:nextPos(1)-1));
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%