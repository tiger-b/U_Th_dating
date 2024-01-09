% ------------------------------------------------------------------- %
% Sub-module that creates an array of the flags for easier export to.
% output. Also prints warning messages if any flags are activated.
% ------------------------------------------------------------------- %
function [flags] = MIT_printFlagMessages(flagStruct)
flags = cell(1, 6); % Pre-allocate space
if ~isequal(flagStruct.Th.varMassBias, 0)
    fprintf('\t\t\t WARNING! Variable Th mass bias. Relative diff = %f per mil.\n', ...
        flagStruct.Th.varMassBias*1000)
    flags(1) = num2cell(flagStruct.Th.varMassBias*1000);
end
if ~isequal(flagStruct.U.varMassBias, 0)
    fprintf('\t\t\t WARNING! Variable U mass bias. Relative diff = %f per mil.\n', ...
        flagStruct.U.varMassBias*1000)
    flags(2) = num2cell(flagStruct.U.varMassBias*1000);
end
if ~isequal(flagStruct.Th.varICYield, 0)
    fprintf('\t\t\t WARNING! Variable Th IC yield. Relative diff = %f per mil.\n', ...
        flagStruct.Th.varICYield*1000)
    flags(3) = num2cell(flagStruct.Th.varICYield*1000);
end
if ~isequal(flagStruct.U.varICYield, 0)
    fprintf('\t\t\t WARNING! Variable U IC yield. Relative diff = %f per mil.\n', ...
        flagStruct.U.varICYield*1000)
    flags(4) = num2cell(flagStruct.U.varICYield*1000);
end
if ~isequal(flagStruct.Th.lowspike, 0)
    fprintf('\t\t\t WARNING! Low Th spike intensity. Relative diff = %f per mil.\n', ...
        flagStruct.Th.lowspike*1000)
    flags(5) = num2cell(flagStruct.Th.lowspike*1000);
end
if ~isequal(flagStruct.U.lowspike, 0)
    fprintf('\t\t\t WARNING! Low U spike intensity. Relative diff = %f per mil.\n', ...
        flagStruct.U.lowspike*1000)
    flags(6) = num2cell(flagStruct.U.lowspike*1000);
end
end