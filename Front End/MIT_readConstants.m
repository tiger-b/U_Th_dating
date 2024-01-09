%% ---------------------- MIT_readConstants.m --------------------- %%
% DESCRIPTION
%   Reads in the values of the following constants:
%       - spike constants (aragonite and calcite)
%       - decay constants and masses
%       - standard values
%
% DEPENDENCIES
%   None
% ----------------------------------------------------------------------- %
% AUTHORS: Christine Y Chen [CYC], David McGee [DM]
% Created:      06/09/15
% Last updated: 06/09/15
%% ---------------------- REVISION HISTORY ------------------------------ %
%
% 06/09/15    [CYC]     Created by CYC. Based on Excel reduction sheet from
%                       Dec 2014.
%
%% --------------------------------------------------------------------- %%
% ---------------------- Read in spike constants ------------------------ %

constants.spk.arag.U238_U233   = 0.000284180083703407;
constants.spk.arag.U235_U233   = 0.0000436564125150588;
constants.spk.arag.U234_U233   = 0.000346078761686624;
constants.spk.arag.Th232_Th229 = 0.0035688899585085;
constants.spk.arag.Th230_Th229 = 0.0000568907919271831;
constants.spk.arag.Th229       = 0.44434;
constants.spk.arag.U233        = 8.25560374325385;
constants.spk.arag.U236_U233   = 0.981296488921163;

constants.spk.calc.U238_U233   = 0.000175450777882931;
constants.spk.calc.U235_U233   = 0.0000430555690434635;
constants.spk.calc.U234_U233   = constants.spk.arag.U234_U233;
constants.spk.calc.Th232_Th229 = 0.00119508371329823;
constants.spk.calc.Th230_Th229 = 0.0000538350424096104;
constants.spk.calc.Th229       = 0.34846;
constants.spk.calc.U233        = 2.572262028;
constants.spk.calc.U236_U233   = 0.981296488921163;

%% --------------------- Read in decay constants and masses ------------- %

constants.decay.U238 = 1.55125E-10;
constants.decay.U234 = 0.00000282206055;
constants.decay.Th232 = 4.9475E-11;
constants.decay.Th230 = 0.00000917052078;
constants.mass.U238 = 238.050785;
constants.mass.U236 = 236.045563;
constants.mass.U235 = 235.043924;
constants.mass.U234 = 234.040947;
constants.mass.U233 = 233.039629;
constants.mass.Th232 = 232.0380553;
constants.mass.Th230 = 230.0331338;
constants.mass.Th229 = 229.031762;

%% --------------------- Read in standard values ------------------------ %

constants.std.CRM112a.U238_U234_mean = 18920.8;
constants.std.CRM112a.U238_U234_sig = 2.7;
constants.std.CRM112a.U238_U235_mean = 137.829;
constants.std.CRM112a.U238_U235_sig = 0.011;
constants.std.IRMM.U236_U233_mean = 0.981296488921163;
constants.std.IRMM.U236_U233_sig = 0.00015;
constants.std.MITh_1.Th232_Th229_mean = 12.3976180730787;
constants.std.MITh_1.Th232_Th229_sig = 0.00203328639774461;
constants.std.MITh_1.Th230_Th229_mean = 0.0503780853190189;
constants.std.MITh_1.Th230_Th229_sig = 4.39868269075332E-06;