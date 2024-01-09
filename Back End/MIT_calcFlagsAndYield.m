%% ------------------------ MIT_calcFlagsAndYield.m --------------------- %%
% DESCRIPTION
%   This module calculates the percent yield of samples and creates flags
%   based on the behavior of the bracketing standards.
%
% INPUTS
%   stanTh1_r:  struct containing all raw Th data for first standard
%   stanTh2_r:  struct containing all raw Th data for second standard
%   sampTh_r:   struct containing all raw Th data for sample 
%   stanU1_r:   struct containing all raw U data for first standard
%   stanU2_r:   struct containing all raw U data for second standard
%   sampU_r:    struct containing all raw U data for sample
%   vol_U:      number of mL of U dissolved fraction
%   vol_Th:     number of mL of Th dissolved fraction
%
% OUTPUTS
%   flag:  struct containing the percent of bracketing standard variation if
%          a specific threshold is passed; otherwise, value is 0
%   yield: struct containing the percent yield in 236U and 229Th
%
% CONSTANTS
%   CONST
% ----------------------------------------------------------------------- %
% AUTHORS: Christine Y Chen [CYC], David McGee [DM]
% Created:      09/04/15
% Last updated: 09/08/15
%% ---------------------- REVISION HISTORY ------------------------------ %
% 90/08/15    [CYC]     Corrected the perent chemistry yield calculations.
% 09/04/15    [CYC]     Created by CYC in order to calculate flags and
%                       percent yield for blank as well as all samples.
%
%% --------------------------------------------------------------------- %%
function [flag, percYield] = MIT_calcFlagsAndYield(stanTh1_r, stanTh2_r, sampTh_r, ...
                                                   stanU1_r,  stanU2_r,  sampU_r, ...
                                                   powder)
% Call global constants
global CONST;

%% ----------------------- Investigate flags --------------------------- %%
% Flag if Th232/Th229 of bracketing standards for a sample change by 
% >0.5 per mil ("variable Th mass bias")
if abs(stanTh1_r.Th232_Th229 - stanTh2_r.Th232_Th229)/stanTh1_r.Th232_Th229 > 0.0005 
    flag.Th.varMassBias = abs(stanTh1_r.Th232_Th229 - stanTh2_r.Th232_Th229)/stanTh1_r.Th232_Th229;
else
    flag.Th.varMassBias = 0;
end

% Flag if U238/U235 of bracketing standards for a sample change by
% >0.5 per mil ("variable U mass bias")
if abs(stanU1_r.U.U238_U235  - stanU2_r.U.U238_U235)/stanU1_r.U.U238_U235 > 0.0005 
    flag.U.varMassBias = abs(stanU1_r.U.U238_U235  - stanU2_r.U.U238_U235)/stanU1_r.U.U238_U235;
else
    flag.U.varMassBias = 0;
end

% Flag if Th230/Th229 of bracketing standards for a sample change by
% >2 per mil ("variable Th IC yield")
if abs(stanTh1_r.Th230_Th229 - stanTh2_r.Th230_Th229)/stanTh1_r.Th230_Th229 > .002
    flag.Th.varICYield = abs(stanTh1_r.Th230_Th229 - stanTh2_r.Th230_Th229)/stanTh1_r.Th230_Th229;
else
    flag.Th.varICYield = 0;
end

% Flag if Th238/Th234 of bracketing standards for a sample change by
% >2 per mil ("variable U IC yield")
if abs(stanU1_r.U.U238_U234  - stanU2_r.U.U238_U234)/stanU1_r.U.U238_U234 > .002 
    flag.U.varICYield = abs(stanU1_r.U.U238_U234  - stanU2_r.U.U238_U234)/stanU1_r.U.U238_U234;
else
    flag.U.varICYield = 0;
end

% Flag if Th229 is <10 mV for a sample ("low Th spike intensity")
if sampTh_r.Th229_V < 0.01 
    flag.Th.lowspike = sampTh_r.Th229_V;
else
    flag.Th.lowspike = 0;
end

% Flag if U236 is <25 mV for a sample ("low U spike intensity")
if sampU_r.U.U236_V < 0.025
    flag.U.lowspike = sampU_r.U.U236_V;
else
    flag.U.lowspike = 0;
end

%% ------------------ Calculate percent yield flags -------------------- %%
% Extract U-Th spike ratios (calcite or aragonite).
if powder.spk_type == 1    % Calcite spike
    spk = CONST.spk.calc;
else                       % Aragonite spike
    spk = CONST.spk.arag;
end

% Calculate percent chemistry yield for U236
U236_expected = powder.spk_wt * spk.U236_U233 * spk.U233 / powder.U_vol * 236 * ...
                   mean([stanU1_r.U.U238_V stanU2_r.U.U238_V])/powder.std_U_conc / 1000; % Sensitivity

percYield(1,1) = sampU_r.U.U236_V/U236_expected*100; 

% Calculate percent chemistry yield for Th229
Th229_expected = powder.spk_wt * spk.Th229 / powder.Th_vol * 229 * ...
                    mean([stanTh1_r.Th232_V stanTh2_r.Th232_V])/powder.std_Th_conc / 1000; % Sensitivity
                    1000; % Convert to volts

percYield(1,2) = sampTh_r.Th229_V/Th229_expected*100;    

end