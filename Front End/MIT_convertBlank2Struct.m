%% ---------------------- MIT_convertBlank2Struct.m --------------------- %%
% DESCRIPTION
%   Converts matrix of blank values into a structure.
%
% CONSTANTS
%   None
% 
% DEPENDENCIES
%   None
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
% ------------------ Put blank values into structure ------------------ %

blank.spk_wt             = blank_spkWt;

% ---------------- Read in Th data for the blank ------------------------ %
blankTh_r.Th230_Th232     = blank_Thsamp_data(1);
blankTh_r.Th230_Th232_err = blank_Thsamp_data(2);
blankTh_r.Th232_Th229     = blank_Thsamp_data(3);
blankTh_r.Th232_Th229_err = blank_Thsamp_data(4);
blankTh_r.Th230_Th229     = blank_Thsamp_data(5);
blankTh_r.Th230_Th229_err = blank_Thsamp_data(6);
blankTh_r.Th232_V         = blank_Thsamp_data(7);
blankTh_r.Th232_V_err     = blank_Thsamp_data(8);
blankTh_r.Th230_cps       = blank_Thsamp_data(9);
blankTh_r.Th230_cps_err   = blank_Thsamp_data(10);
blankTh_r.Th229_V         = blank_Thsamp_data(11);
blankTh_r.Th229_V_err     = blank_Thsamp_data(12);
blankTh_r.U238_V          = blank_Thsamp_data(13);
blankTh_r.U238_V_err      = blank_Thsamp_data(14);
blankTh_r.Th230_bg        = blank_Thsamp_data(15);
blankTh_r.Th230_bg_err    = blank_Thsamp_data(16);

% ----------- Read in Th data for Standard #1 of blank ------------------ %
stdBlnkTh1_r.Th230_Th232     = blank_Thstd1_data(1);
stdBlnkTh1_r.Th230_Th232_err = blank_Thstd1_data(2);
stdBlnkTh1_r.Th232_Th229     = blank_Thstd1_data(3);
stdBlnkTh1_r.Th232_Th229_err = blank_Thstd1_data(4);
stdBlnkTh1_r.Th230_Th229     = blank_Thstd1_data(5);
stdBlnkTh1_r.Th230_Th229_err = blank_Thstd1_data(6);
stdBlnkTh1_r.Th232_V         = blank_Thstd1_data(7);
stdBlnkTh1_r.Th232_V_err     = blank_Thstd1_data(8);
stdBlnkTh1_r.Th230_cps       = blank_Thstd1_data(9);
stdBlnkTh1_r.Th230_cps_err   = blank_Thstd1_data(10);
stdBlnkTh1_r.Th229_V         = blank_Thstd1_data(11);
stdBlnkTh1_r.Th229_V_err     = blank_Thstd1_data(12);
stdBlnkTh1_r.U238_V          = blank_Thstd1_data(13);
stdBlnkTh1_r.U238_V_err      = blank_Thstd1_data(14);
stdBlnkTh1_r.Th230_bg        = blank_Thstd1_data(15);
stdBlnkTh1_r.Th230_bg_err    = blank_Thstd1_data(16);

% ------------ Read in Th data for Standard #2 of blank------------------ %
stdBlnkTh2_r.Th230_Th232     = blank_Thstd2_data(1);
stdBlnkTh2_r.Th230_Th232_err = blank_Thstd2_data(2);
stdBlnkTh2_r.Th232_Th229     = blank_Thstd2_data(3);
stdBlnkTh2_r.Th232_Th229_err = blank_Thstd2_data(4);
stdBlnkTh2_r.Th230_Th229     = blank_Thstd2_data(5);
stdBlnkTh2_r.Th230_Th229_err = blank_Thstd2_data(6);
stdBlnkTh2_r.Th232_V         = blank_Thstd2_data(7);
stdBlnkTh2_r.Th232_V_err     = blank_Thstd2_data(8);
stdBlnkTh2_r.Th230_cps       = blank_Thstd2_data(9);
stdBlnkTh2_r.Th230_cps_err   = blank_Thstd2_data(10);
stdBlnkTh2_r.Th229_V         = blank_Thstd2_data(11);
stdBlnkTh2_r.Th229_V_err     = blank_Thstd2_data(12);
stdBlnkTh2_r.U238_V          = blank_Thstd2_data(13);
stdBlnkTh2_r.U238_V_err      = blank_Thstd2_data(14);
stdBlnkTh2_r.Th230_bg        = blank_Thstd2_data(15);
stdBlnkTh2_r.Th230_bg_err    = blank_Thstd2_data(16);

% ---------------- Read in Standard 1 U of the blank -------------------- %
stdBlnkU1_r.U.U238_U236      = blank_Ustd1_data(1);
stdBlnkU1_r.U.U238_U236_err  = blank_Ustd1_data(2);
stdBlnkU1_r.U.U238_U235      = blank_Ustd1_data(3);
stdBlnkU1_r.U.U238_U235_err  = blank_Ustd1_data(4);
stdBlnkU1_r.U.U236_U233      = blank_Ustd1_data(5);
stdBlnkU1_r.U.U236_U233_err  = blank_Ustd1_data(6);
stdBlnkU1_r.U.U238_U234      = blank_Ustd1_data(7);
stdBlnkU1_r.U.U238_U234_err  = blank_Ustd1_data(8);
stdBlnkU1_r.U.U238_V         = blank_Ustd1_data(9);
stdBlnkU1_r.U.U238_V_err     = blank_Ustd1_data(10);
stdBlnkU1_r.U.U236_V         = blank_Ustd1_data(11);
stdBlnkU1_r.U.U236_V_err     = blank_Ustd1_data(12);
stdBlnkU1_r.U.U234_cps       = blank_Ustd1_data(13);
stdBlnkU1_r.U.U234_cps_err   = blank_Ustd1_data(14);

% -------------- Read in Standard 1 U tail of the blank ----------------- %
stdBlnkU1_r.Utail.U238_V       = blank_Ustd1_tail_data(1);
stdBlnkU1_r.Utail.U238_V_err   = blank_Ustd1_tail_data(2);
stdBlnkU1_r.Utail.U236_V       = blank_Ustd1_tail_data(3);
stdBlnkU1_r.Utail.U236_V_err   = blank_Ustd1_tail_data(4);
stdBlnkU1_r.Utail.U237_cps     = blank_Ustd1_tail_data(5);
stdBlnkU1_r.Utail.U237_err     = blank_Ustd1_tail_data(6);
stdBlnkU1_r.Utail.U2365_cps    = blank_Ustd1_tail_data(7);
stdBlnkU1_r.Utail.U2365_err    = blank_Ustd1_tail_data(8);
stdBlnkU1_r.Utail.U2345_cps    = blank_Ustd1_tail_data(9);
stdBlnkU1_r.Utail.U2345_err    = blank_Ustd1_tail_data(10);
stdBlnkU1_r.Utail.U2335_cps    = blank_Ustd1_tail_data(11);
stdBlnkU1_r.Utail.U2335_err    = blank_Ustd1_tail_data(12);

% ---------------- Read in Standard 2 U tail of the blank --------------- %
stdBlnkU2_r.Utail.U238_V       = blank_Ustd2_tail_data(1);
stdBlnkU2_r.Utail.U238_V_err   = blank_Ustd2_tail_data(2);
stdBlnkU2_r.Utail.U236_V       = blank_Ustd2_tail_data(3);
stdBlnkU2_r.Utail.U236_V_err   = blank_Ustd2_tail_data(4);
stdBlnkU2_r.Utail.U237_cps     = blank_Ustd2_tail_data(5);
stdBlnkU2_r.Utail.U237_err     = blank_Ustd2_tail_data(6);
stdBlnkU2_r.Utail.U2365_cps    = blank_Ustd2_tail_data(7);
stdBlnkU2_r.Utail.U2365_err    = blank_Ustd2_tail_data(8);
stdBlnkU2_r.Utail.U2345_cps    = blank_Ustd2_tail_data(9);
stdBlnkU2_r.Utail.U2345_err    = blank_Ustd2_tail_data(10);
stdBlnkU2_r.Utail.U2335_cps    = blank_Ustd2_tail_data(11);
stdBlnkU2_r.Utail.U2335_err    = blank_Ustd2_tail_data(12);

% --------------------- Read in U data for the blank -------------------- %
blankU_r.U.U238_U236      = blank_Usamp_data(1);
blankU_r.U.U238_U236_err  = blank_Usamp_data(2);
blankU_r.U.U238_U235      = blank_Usamp_data(3);
blankU_r.U.U238_U235_err  = blank_Usamp_data(4);
blankU_r.U.U236_U233      = blank_Usamp_data(5);
blankU_r.U.U236_U233_err  = blank_Usamp_data(6);
blankU_r.U.U238_U234      = blank_Usamp_data(7);
blankU_r.U.U238_U234_err  = blank_Usamp_data(8);
blankU_r.U.U238_V         = blank_Usamp_data(9);
blankU_r.U.U238_V_err     = blank_Usamp_data(10);
blankU_r.U.U236_V         = blank_Usamp_data(11);
blankU_r.U.U236_V_err     = blank_Usamp_data(12);
blankU_r.U.U234_cps       = blank_Usamp_data(13);
blankU_r.U.U234_cps_err   = blank_Usamp_data(14);

% -------------------- Read in U tail of the blank ---------------------- %
blankU_r.Utail.U238_V        = blank_Usamp_data(9);
blankU_r.Utail.U238_V_err    = blank_Usamp_data(10);
blankU_r.Utail.U236_V        = blank_Usamp_data(11);
blankU_r.Utail.U236_V_err    = blank_Usamp_data(12);
%blankU_r.Utail.U237_cps      = blank_Usamp_tail_data(5);
blankU_r.Utail.U237_cps      = ((blank_Ustd1_tail_data(5)/blank_Ustd1_tail_data(1)) + (blank_Ustd2_tail_data(5)/blank_Ustd2_tail_data(1)))/2 * blank_Usamp_data(9);
%blankU_r.Utail.U237_err      = blank_Usamp_tail_data(6);
blankU_r.Utail.U237_err       = ((blank_Ustd1_tail_data(6)/blank_Ustd1_tail_data(5)) + (blank_Ustd2_tail_data(6)/blank_Ustd2_tail_data(5)))/2 * blankU_r.Utail.U237_cps;
%blankU_r.Utail.U2365_cps     = blank_Usamp_tail_data(7);
blankU_r.Utail.U2365_cps     = ((blank_Ustd1_tail_data(7)/blank_Ustd1_tail_data(1)) + (blank_Ustd2_tail_data(7)/blank_Ustd2_tail_data(1)))/2 * blank_Usamp_data(9);
%blankU_r.Utail.U2365_err     = blank_Usamp_tail_data(8);
blankU_r.Utail.U2365_err     = ((blank_Ustd1_tail_data(8)/blank_Ustd1_tail_data(7)) + (blank_Ustd2_tail_data(8)/blank_Ustd2_tail_data(7)))/2 * blankU_r.Utail.U2365_cps;
%blankU_r.Utail.U2345_cps     = blank_Usamp_tail_data(9);
blankU_r.Utail.U2345_cps     = ((blank_Ustd1_tail_data(9)/blank_Ustd1_tail_data(1)) + (blank_Ustd2_tail_data(9)/blank_Ustd2_tail_data(1)))/2 * blank_Usamp_data(9);
%blankU_r.Utail.U2345_err     = blank_Usamp_tail_data(10);
blankU_r.Utail.U2345_err     = ((blank_Ustd1_tail_data(10)/blank_Ustd1_tail_data(9)) + (blank_Ustd2_tail_data(10)/blank_Ustd2_tail_data(9)))/2 * blankU_r.Utail.U2345_cps;
%blankU_r.Utail.U2335_cps     = blank_Usamp_tail_data(11);
blankU_r.Utail.U2335_cps     = ((blank_Ustd1_tail_data(11)/blank_Ustd1_tail_data(1)) + (blank_Ustd2_tail_data(11)/blank_Ustd2_tail_data(1)))/2 * blank_Usamp_data(9);
%blankU_r.Utail.U2335_err     = blank_Usamp_tail_data(12);
blankU_r.Utail.U2335_err     = ((blank_Ustd1_tail_data(12)/blank_Ustd1_tail_data(11)) + (blank_Ustd2_tail_data(12)/blank_Ustd2_tail_data(11)))/2 * blankU_r.Utail.U2335_cps;


% ----------------- Read in Standard 2 U of the blank-------------------- %
stdBlnkU2_r.U.U238_U236      = blank_Ustd2_data(1);
stdBlnkU2_r.U.U238_U236_err  = blank_Ustd2_data(2);
stdBlnkU2_r.U.U238_U235      = blank_Ustd2_data(3);
stdBlnkU2_r.U.U238_U235_err  = blank_Ustd2_data(4);
stdBlnkU2_r.U.U236_U233      = blank_Ustd2_data(5);
stdBlnkU2_r.U.U236_U233_err  = blank_Ustd2_data(6);
stdBlnkU2_r.U.U238_U234      = blank_Ustd2_data(7);
stdBlnkU2_r.U.U238_U234_err  = blank_Ustd2_data(8);
stdBlnkU2_r.U.U238_V         = blank_Ustd2_data(9);
stdBlnkU2_r.U.U238_V_err     = blank_Ustd2_data(10);
stdBlnkU2_r.U.U236_V         = blank_Ustd2_data(11);
stdBlnkU2_r.U.U236_V_err     = blank_Ustd2_data(12);
stdBlnkU2_r.U.U234_cps       = blank_Ustd2_data(13);
stdBlnkU2_r.U.U234_cps_err   = blank_Ustd2_data(14);



% Clear blank data that is not in a structure
clear blank_Thsamp_data blank_Thstd1_data blank_Thstd2_data blank_Usamp_data ...     
      blank_Usamp_tail_data blank_Ustd1_data blank_Ustd1_tail_data blank_Ustd2_data blank_Ustd2_tail_data 
