%% ---------------------- MIT_convertSamp2Struct.m --------------------- %%
% DESCRIPTION
%   Converts sample matrix data into a structure.
%
% INPUTS
%   i      : The number of the sample being looked at from loop
%   Thsamp : Th sample data matrix (i number of rows)
%   Thstd  : Th standard data matrix
%   Usamp  : U sample data matrix
%   Ustd   : U standard data matrix
%   Utail  : U tail data matrix
% OUTPUTS
%   sampleTh    : structure of one sample Th data
%   standardTh1 : structure of standard Th 1 data
%   standardTh2 : structure of standard Th 2 data
%   sampleU     : structure of one sample U data
%   standardU1  : structure of standard U1 data
%   standardU2  : strucutre of standard U2 data
%
% CONSTANTS
%   None.
% 
% DEPENDENCIES
%   None.
% ----------------------------------------------------------------------- %
% AUTHORS: Christine Y Chen [CYC], David McGee [DM]
% Created:      06/12/15
% Last updated: 06/12/15
%% ---------------------- REVISION HISTORY ------------------------------ %
% 08/13/15    [CYC]     Changes reflected due to the new raw summary file 
%                       for U-tail having fewer columns.
% 06/12/15    [CYC]     Created by CYC. Based on Excel reduction sheet from
%                       Dec 2014.
%
%% --------------------------------------------------------------------- %%
function [sampleTh, standardTh1, standardTh2, sampleU, standardU1, standardU2] ...
        = MIT_convertSamp2Struct(i, Thsamp, Thstd1, Thstd2, ...
                                    Usamp, Usamp_tail, Ustd1, Ustd1_tail, ...
                                    Ustd2, Ustd2_tail)

    % ---------------- Read in Th data for sample ----------------------- %
    sampleTh.Th230_Th232     = Thsamp(i,  1);
    sampleTh.Th230_Th232_err = Thsamp(i,  2);
    sampleTh.Th232_Th229     = Thsamp(i,  3);
    sampleTh.Th232_Th229_err = Thsamp(i,  4);
    sampleTh.Th230_Th229     = Thsamp(i,  5);
    sampleTh.Th230_Th229_err = Thsamp(i,  6);
    sampleTh.Th232_V         = Thsamp(i,  7);
    sampleTh.Th232_V_err     = Thsamp(i,  8);
    sampleTh.Th230_cps       = Thsamp(i,  9);
    sampleTh.Th230_cps_err   = Thsamp(i, 10);
    sampleTh.Th229_V         = Thsamp(i, 11);
    sampleTh.Th229_V_err     = Thsamp(i, 12);
    sampleTh.U238_V          = Thsamp(i, 13);
    sampleTh.U238_V_err      = Thsamp(i, 14);
    sampleTh.Th230_bg        = Thsamp(i, 15);
    sampleTh.Th230_bg_err    = Thsamp(i, 16);
    
    % ---------------- Read in Th data for Standard #1 ------------------ %
    standardTh1.Th230_Th232     = Thstd1(i,  1);
    standardTh1.Th230_Th232_err = Thstd1(i,  2);
    standardTh1.Th232_Th229     = Thstd1(i,  3);
    standardTh1.Th232_Th229_err = Thstd1(i,  4);
    standardTh1.Th230_Th229     = Thstd1(i,  5);
    standardTh1.Th230_Th229_err = Thstd1(i,  6);
    standardTh1.Th232_V         = Thstd1(i,  7);
    standardTh1.Th232_V_err     = Thstd1(i,  8);
    standardTh1.Th230_cps       = Thstd1(i,  9);
    standardTh1.Th230_cps_err   = Thstd1(i, 10);
    standardTh1.Th229_V         = Thstd1(i, 11);
    standardTh1.Th229_V_err     = Thstd1(i, 12);
    standardTh1.U238_V          = Thstd1(i, 13);
    standardTh1.U238_V_err      = Thstd1(i, 14);
    standardTh1.Th230_bg        = Thstd1(i, 15);
    standardTh1.Th230_bg_err    = Thstd1(i, 16);
     
    % ---------------- Read in Th data for Standard #2 ------------------ %
    standardTh2.Th230_Th232     = Thstd2(i,  1);
    standardTh2.Th230_Th232_err = Thstd2(i,  2);
    standardTh2.Th232_Th229     = Thstd2(i,  3);
    standardTh2.Th232_Th229_err = Thstd2(i,  4);
    standardTh2.Th230_Th229     = Thstd2(i,  5);
    standardTh2.Th230_Th229_err = Thstd2(i,  6);
    standardTh2.Th232_V         = Thstd2(i,  7);
    standardTh2.Th232_V_err     = Thstd2(i,  8);
    standardTh2.Th230_cps       = Thstd2(i,  9);
    standardTh2.Th230_cps_err   = Thstd2(i, 10);
    standardTh2.Th229_V         = Thstd2(i, 11);
    standardTh2.Th229_V_err     = Thstd2(i, 12);
    standardTh2.U238_V          = Thstd2(i, 13);
    standardTh2.U238_V_err      = Thstd2(i, 14);
    standardTh2.Th230_bg        = Thstd2(i, 15);
    standardTh2.Th230_bg_err    = Thstd2(i, 16);
    
    % ---------------- Standard 1 U ------------------------------------- %
    standardU1.U.U238_U236      = Ustd1(i,  1);
    standardU1.U.U238_U236_err  = Ustd1(i,  2);
    standardU1.U.U238_U235      = Ustd1(i,  3);
    standardU1.U.U238_U235_err  = Ustd1(i,  4);
    standardU1.U.U236_U233      = Ustd1(i,  5);
    standardU1.U.U236_U233_err  = Ustd1(i,  6);
    standardU1.U.U238_U234      = Ustd1(i,  7);
    standardU1.U.U238_U234_err  = Ustd1(i,  8);
    standardU1.U.U238_V         = Ustd1(i,  9);
    standardU1.U.U238_V_err     = Ustd1(i, 10);
    standardU1.U.U236_V         = Ustd1(i, 11);
    standardU1.U.U236_V_err     = Ustd1(i, 12);
    standardU1.U.U234_cps       = Ustd1(i, 13);
    standardU1.U.U234_cps_err   = Ustd1(i, 14);
    
    % ------------------ Standard 1 U tail ------------------------------ %
    standardU1.Utail.U238_V       = Ustd1_tail(i, 1);
    standardU1.Utail.U238_V_err   = Ustd1_tail(i, 2);
    standardU1.Utail.U236_V       = Ustd1_tail(i, 3);
    standardU1.Utail.U236_V_err   = Ustd1_tail(i, 4);
    standardU1.Utail.U237_cps     = Ustd1_tail(i, 5);
    standardU1.Utail.U237_err     = Ustd1_tail(i, 6);
    standardU1.Utail.U2365_cps    = Ustd1_tail(i, 7);
    standardU1.Utail.U2365_err    = Ustd1_tail(i, 8);
    standardU1.Utail.U2345_cps    = Ustd1_tail(i, 9);
    standardU1.Utail.U2345_err    = Ustd1_tail(i, 10);
    standardU1.Utail.U2335_cps    = Ustd1_tail(i, 11);
    standardU1.Utail.U2335_err    = Ustd1_tail(i, 12);
    
    % ------------------ Standard 2 U tail ------------------------------ %
    standardU2.Utail.U238_V       = Ustd2_tail(i, 1);
    standardU2.Utail.U238_V_err   = Ustd2_tail(i, 2);
    standardU2.Utail.U236_V       = Ustd2_tail(i, 3);
    standardU2.Utail.U236_V_err   = Ustd2_tail(i, 4);
    standardU2.Utail.U237_cps     = Ustd2_tail(i, 5);
    standardU2.Utail.U237_err     = Ustd2_tail(i, 6);
    standardU2.Utail.U2365_cps    = Ustd2_tail(i, 7);
    standardU2.Utail.U2365_err    = Ustd2_tail(i, 8);
    standardU2.Utail.U2345_cps    = Ustd2_tail(i, 9);
    standardU2.Utail.U2345_err    = Ustd2_tail(i, 10);
    standardU2.Utail.U2335_cps    = Ustd2_tail(i, 11);
    standardU2.Utail.U2335_err    = Ustd2_tail(i, 12);

    % -------------------- Read in U data for sample -------------------- %
    sampleU.U.U238_U236      = Usamp(i,  1);
    sampleU.U.U238_U236_err  = Usamp(i,  2);
    sampleU.U.U238_U235      = Usamp(i,  3);
    sampleU.U.U238_U235_err  = Usamp(i,  4);
    sampleU.U.U236_U233      = Usamp(i,  5);
    sampleU.U.U236_U233_err  = Usamp(i,  6);
    sampleU.U.U238_U234      = Usamp(i,  7);
    sampleU.U.U238_U234_err  = Usamp(i,  8);
    sampleU.U.U238_V         = Usamp(i,  9);
    sampleU.U.U238_V_err     = Usamp(i, 10);
    sampleU.U.U236_V         = Usamp(i, 11);
    sampleU.U.U236_V_err     = Usamp(i, 12);
    sampleU.U.U234_cps       = Usamp(i, 13);
    sampleU.U.U234_cps_err   = Usamp(i, 14);
    
    % ------------------------ Sample U tail ---------------------------- %
    %sampleU.Utail.U238_V        = Usamp_tail(i,  1);
    sampleU.Utail.U238_V         = Usamp(i,  9);
    %sampleU.Utail.U238_V_err    = Usamp_tail(i,  2);
    sampleU.Utail.U238_V_err    = Usamp(i, 10);
    %sampleU.Utail.U236_V        = Usamp_tail(i,  3);
    sampleU.Utail.U236_V        = Usamp(i, 11);
    %sampleU.Utail.U236_V_err    = Usamp_tail(i,  4);
    sampleU.Utail.U236_V_err    = Usamp(i, 12);
    %sampleU.Utail.U237_cps      = Usamp_tail(i,  5);
    sampleU.Utail.U237_cps      = ((Ustd1_tail(i, 5)/Ustd1_tail(i, 1)) + (Ustd2_tail(i, 5)/Ustd2_tail(i, 1)))/2 * Usamp(i,  9);
    %sampleU.Utail.U237_err      = Usamp_tail(i,  6);
    sampleU.Utail.U237_err      = ((Ustd1_tail(i, 6)/Ustd1_tail(i, 5)) + (Ustd2_tail(i, 6)/Ustd2_tail(i, 5)))/2 * sampleU.Utail.U237_cps;
    %sampleU.Utail.U2365_cps     = Usamp_tail(i,  7);
    sampleU.Utail.U2365_cps     = ((Ustd1_tail(i, 7)/Ustd1_tail(i, 1)) + (Ustd2_tail(i, 7)/Ustd2_tail(i, 1)))/2 * Usamp(i,  9);
    %sampleU.Utail.U2365_err     = Usamp_tail(i,  8);
    sampleU.Utail.U2365_err     = ((Ustd1_tail(i, 8)/Ustd1_tail(i, 7)) + (Ustd2_tail(i, 8)/Ustd2_tail(i, 7)))/2 * sampleU.Utail.U2365_cps;
    %sampleU.Utail.U2345_cps     = Usamp_tail(i,  9);
    sampleU.Utail.U2345_cps     = ((Ustd1_tail(i, 9)/Ustd1_tail(i, 1)) + (Ustd2_tail(i, 9)/Ustd2_tail(i, 1)))/2 * Usamp(i,  9);
    %sampleU.Utail.U2345_err     = Usamp_tail(i, 10);
    sampleU.Utail.U2345_err     = ((Ustd1_tail(i, 10)/Ustd1_tail(i, 9)) + (Ustd2_tail(i, 10)/Ustd2_tail(i, 9)))/2 * sampleU.Utail.U2345_cps;
    %sampleU.Utail.U2335_cps     = Usamp_tail(i, 11);
    sampleU.Utail.U2335_cps     = ((Ustd1_tail(i, 11)/Ustd1_tail(i, 1)) + (Ustd2_tail(i, 11)/Ustd2_tail(i, 1)))/2 * Usamp(i,  9);
    %sampleU.Utail.U2335_err     = Usamp_tail(i, 12);
    sampleU.Utail.U2335_err     = ((Ustd1_tail(i, 12)/Ustd1_tail(i, 11)) + (Ustd2_tail(i, 12)/Ustd2_tail(i, 11)))/2 * sampleU.Utail.U2335_cps;
    
    % ---------------- Standard 2 U ------------------------------------- %
    standardU2.U.U238_U236      = Ustd2(i,  1);
    standardU2.U.U238_U236_err  = Ustd2(i,  2);
    standardU2.U.U238_U235      = Ustd2(i,  3);
    standardU2.U.U238_U235_err  = Ustd2(i,  4);
    standardU2.U.U236_U233      = Ustd2(i,  5);
    standardU2.U.U236_U233_err  = Ustd2(i,  6);
    standardU2.U.U238_U234      = Ustd2(i,  7);
    standardU2.U.U238_U234_err  = Ustd2(i,  8);
    standardU2.U.U238_V         = Ustd2(i,  9);
    standardU2.U.U238_V_err     = Ustd2(i, 10);
    standardU2.U.U236_V         = Ustd2(i, 11);
    standardU2.U.U236_V_err     = Ustd2(i, 12);
    standardU2.U.U234_cps       = Ustd2(i, 13);
    standardU2.U.U234_cps_err   = Ustd2(i, 14);
    

end