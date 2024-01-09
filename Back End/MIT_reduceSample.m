%% ---------------------- MIT_reduceSample.m --------------------- %%
% DESCRIPTION
%   Workhorse of the U-Th data reduction scheme. Corrects one sample
%   bracketed by two standards.
%
% INPUTS
%   stanTh1_r:  struct containing all raw Th data for first standard
%   stanTh2_r:  struct containing all raw Th data for second standard
%   sampTh_r:   struct containing all raw Th data for sample 
%   stanU1_r:   struct containing all raw U data for first standard
%   stanU2_r:   struct containing all raw U data for second standard
%   sampU_r:    struct containing all raw U data for sample 
%   powder:     struct containing info about the sample powder
%   blank_r:    struct containing info about the procedural blank
%
% OUTPUTS
%   Fsamp:      struct of U-Th concentrations for reporting
%   Fages:      struct of final U-Th ages for reporting
%   blank_proc: struct of procedural blank information, including masses
%               and its impact on the final U and Th concentrations
%   flag:       struct of flags for bracketing standards
%   U234t_234:  value of 234tail/234 in per mil
%
% CONSTANTS
%   CONST:      struct of constants used for analysis (e.g., decay
%                  constants, spike U/Th ratios)
% 
% DEPENDENCIES
%   MIT_correctTailMassBiasYield : Tail, mass bias, and yield corrections  
%       MIT_Th_tail_corr.m       : Th tail correction
%       MIT_U_tail_corr.m        : Part 1 of U tail correction
%       MIT_U_tail_corr_beta.m   : Part 2 of U tail correction
%   MIT_spike_corr.m             : Spike correction
%   MIT_Age_calc.m               : Final age calculation
% ----------------------------------------------------------------------- %
% AUTHORS: Christine Y Chen [CYC], David McGee [DM]
% Created:      04/30/15
% Last updated: 09/04/15
%% ---------------------- REVISION HISTORY ------------------------------ %
% 09/04/15    [CYC]     Removed flags and put them in a separate function.
% 08/31/15    [CYC]     Changed procedural blank correction to be uniformly
%                       25% 1-sigma.
%                       Removed 0.95 factor in spike correction of
%                       procedural blank. Now blanks and samples are
%                       spike corrected the same way.
% 08/08/15    [CYC]     Changed flags to output the relative difference of
%                       the two bracketing standards. Changed flags to be
%                       Th and U specific.
% 08/05/15    [CYC]     Added field to procedural blank about its impact on
%                       the final U and Th concentrations.
% 06/12/15    [CYC]     Added section for flags.
% 06/12/15    [CYC]     Changed mass weights in procedural blank correction
%                       (and only in that section) back to integer masses 
%                       in order to be consistent with spreadsheet (for now). 
% 06/09/15    [CYC]     Changed function name to MIT_UTh_reduce_sample.m
%                       instead of MIT_UTh_data_reduction.m.
% 05/01/15    [CYC]     Changed all mass weights (e.g. 234, 238, etc.) in
%                       equations of procedural blank correction to equal 
%                       the true atomic masses as set by the constants.
% 05/01/15    [CYC]     Used global variable CONST to relay constants
%                       instead of having the function take the const 
%                       structure as an input.
% 04/29/15    [CYC]     Created by CYC. Based on Excel reduction sheet from
%                       Dec 2014.
%% --------------------------------------------------------------------- %%

function [Fsamp, Fages, blank_mass, blank_impact, U234t_234] = ...
                      MIT_reduceSample(stanTh1_r, stanTh2_r, sampTh_r, ...
                                       stanU1_r,  stanU2_r,  sampU_r, ...
                                       powder, blank_r)
% Call global constants
global CONST;

%% ------- Tail, mass bias, and yield corrections for U and Th --------- %%

% Output is U and Th concentrations with tail, mass biased, and yield
% corrections (tmby)
[sampTh_tmby, sampU_tmby, U234t_234] = MIT_correctTailMassBiasYield(stanTh1_r, stanTh2_r, sampTh_r, ...
                                                                    stanU1_r,  stanU2_r,  sampU_r);

%% ------------------------ Spike correction --------------------------- %%
% Extract U-Th spike ratios (calcite or aragonite).
if powder.spk_type == 1    % Calcite spike
    spk = CONST.spk.calc;
else                       % Aragonite spike
    spk = CONST.spk.arag;
end

% Correct the blank 
blankOrNot = 1;
[blank_spk.U, blank_spk.Th] = MIT_spike_corr(blank_r.U, blank_r.Th, spk, blankOrNot);
clear blankOrNot

% Then, correct the sample
blankOrNot = 0;
[sampU_spk, sampTh_spk] = MIT_spike_corr(sampU_tmby, sampTh_tmby, spk, blankOrNot);

    % ------------------------------------------------------------------- %
    % Sub-module that applies the spike correction to data.
    % ------------------------------------------------------------------- % 
    function [cor_U, cor_Th] = MIT_spike_corr(raw_U, raw_Th, spike, blankOrNot)
        % Calculate spike correction factors
        cfactor.U238  = 1 - spike.U238_U233/CONST.std.IRMM.U236_U233_mean/raw_U.U238_U236;
        cfactor.Th232 = 1 - spike.Th232_Th229/raw_Th.Th232_Th229;
        cfactor.Th230 = 1 - spike.Th230_Th229/raw_Th.Th230_Th229;

        if blankOrNot % Spike correcting a blank
            cfactor.U234  = 1 - spike.U234_U233*raw_U.U238_U234/raw_U.U238_U236/spike.U236_U233;
            %cfactor.U234  = 1 - spike.U234_U233*(raw_U.U238_U234*0.95)/raw_U.U238_U236/spike.U236_U233;
        else          % Spike correcting a sample
            cfactor.U234  = 1 - spike.U234_U233*raw_U.U238_U234/raw_U.U238_U236/spike.U236_U233;
        end

        % Calculate spike-corrected ratios
        cor_U.U238_U236       = cfactor.U238*raw_U.U238_U236;
        cor_U.U238_U236_err   = sqrt(raw_U.U238_U236_err^2 + ...
                                    (raw_U.U238_U236*(0.001*(1 - cfactor.U238)))^2);
        cor_U.U238_U234       = raw_U.U238_U234*cfactor.U238/cfactor.U234;
        cor_U.U238_U234_err   = sqrt((raw_U.U238_U234_err*cfactor.U238/cfactor.U234)^2 + ...
                                    (raw_U.U238_U234*cfactor.U238*(1 - cfactor.U234)* ...
                                    0.001/cfactor.U234^2)^2);
        cor_Th.Th232_Th229     = cfactor.Th232*raw_Th.Th232_Th229;
        cor_Th.Th232_Th229_err = sqrt((raw_Th.Th232_Th229_err*cfactor.Th232)^2 + ...
                                     (raw_Th.Th232_Th229*(1 - cfactor.Th232)*0.01)^2);
        cor_Th.Th230_Th229     = cfactor.Th230*raw_Th.Th230_Th229;
        cor_Th.Th230_Th229_err = sqrt((raw_Th.Th230_Th229_err*cfactor.Th230)^2 + ...
                                     (raw_Th.Th230_Th229*(1 - cfactor.Th230)*0.005)^2);
        cor_U.U234_U233        = cor_U.U238_U236/cor_U.U238_U234/spike.U236_U233;
    end

%% ------------------ Procedural blank correction ---------------------- %%

% Calculate masses for procedural blank
blank_mass.Th230_fg = blank_spk.Th.Th230_Th229 * spk.Th229 * 230 * blank_r.spk_wt * 1000;
blank_mass.Th232_pg = blank_spk.Th.Th232_Th229 * spk.Th229 * 232 * blank_r.spk_wt;
blank_mass.U234_fg  = blank_spk.U.U234_U233    * spk.U233  * 234 * blank_r.spk_wt * 1000;
blank_mass.U238_pg  = blank_spk.U.U238_U236    * spk.U233  * 238 * blank_r.spk_wt * spk.U236_U233;

% Apply blank correction to sample
sampU_blnk.U238_U236          = sampU_spk.U238_U236*(1 - blank_mass.U238_pg /...
                                (spk.U233 * powder.spk_wt * 238 * sampU_spk.U238_U236 / spk.U236_U233));
sampU_blnk.U238_U236_err      = sqrt(sampU_spk.U238_U236_err^2 + ...
                                    ((sampU_spk.U238_U236 - sampU_blnk.U238_U236)/4)^2);
sampU_blnk.U238_U234          = sampU_spk.U238_U234 * ...
                                (sampU_blnk.U238_U236/sampU_spk.U238_U236)/...
                                (1 - blank_mass.U234_fg/...
                                (spk.U233 * sampU_spk.U234_U233 * powder.spk_wt * 1000 * 234));
sampU_blnk.U238_U234_err      = sqrt(sampU_spk.U238_U234_err^2 + ...
                                    ((sampU_blnk.U238_U234 - sampU_spk.U238_U234)/4)^2); 
sampTh_blnk.Th232_Th229       = sampTh_spk.Th232_Th229  * (1 - blank_mass.Th232_pg/ ...
                                (sampTh_spk.Th232_Th229 * powder.spk_wt * spk.Th229 * 229)); 
sampTh_blnk.Th232_Th229_err   = sqrt(sampTh_spk.Th232_Th229_err^2 + ...
                                ((sampTh_blnk.Th232_Th229 - sampTh_spk.Th232_Th229)/4)^2);
sampTh_blnk.Th230_Th229       = sampTh_spk.Th230_Th229  * (1 - blank_mass.Th230_fg/ ...
                                (sampTh_spk.Th230_Th229 * powder.spk_wt * spk.Th229 * 230 * 1000));
sampTh_blnk.Th230_Th229_err   = sqrt(sampTh_spk.Th230_Th229_err^2 + ...
                                ((sampTh_blnk.Th230_Th229 - sampTh_spk.Th230_Th229)/4)^2);
sampU_blnk.d234U              = (1/sampU_blnk.U238_U234 *(CONST.decay.U234/CONST.decay.U238) - 1) * 1000;
sampU_blnk.d234U_err          = sampU_blnk.U238_U234_err/sampU_blnk.U238_U234 * 1000;
sampU_blnk.Th230_U238_atm     = sampTh_blnk.Th230_Th229/(sampU_blnk.U238_U236 * spk.U236_U233) * spk.Th229/spk.U233;
sampU_blnk.Th230_U238_atm_err = sampU_blnk.Th230_U238_atm * ...
                                sqrt((sampTh_blnk.Th230_Th229_err/sampTh_blnk.Th230_Th229)^2 + ...
                                (sampU_blnk.U238_U236_err/sampU_blnk.U238_U236)^2 + (1/1000)^2);

%% ---------- Final sample concentrations for reporting ---------------- %%
Fsamp.U238             = sampU_blnk.U238_U236 * spk.U233 / spk.U236_U233 * ...
                         powder.spk_wt/powder.samp_wt * CONST.mass.U238 / 1000;  % ng/g
Fsamp.U238_2sig        = 2 * Fsamp.U238 * ...
                         sqrt((sampU_blnk.U238_U236_err/1000)^2 + (0.01)^2);  
Fsamp.Th232            = sampTh_blnk.Th232_Th229 * spk.Th229 * CONST.mass.Th232 * ...
                         powder.spk_wt/powder.samp_wt;  % pg/g
Fsamp.Th232_2sig       = 2 * Fsamp.Th232 * ...
                         sqrt((sampTh_blnk.Th232_Th229_err/sampTh_blnk.Th232_Th229)^2 + 0.01^2);
Fsamp.d234U            = sampU_blnk.d234U;  % per mille
Fsamp.d234U_2sig       = 2 * sampU_blnk.d234U_err;
Fsamp.Th230_U238       = sampU_blnk.Th230_U238_atm * CONST.decay.Th230/CONST.decay.U238;  % activity
Fsamp.Th230_U238_2sig  = 2 * Fsamp.Th230_U238 * sampU_blnk.Th230_U238_atm_err/sampU_blnk.Th230_U238_atm;  
Fsamp.Th230_Th232      = sampTh_blnk.Th230_Th229/sampTh_blnk.Th232_Th229*1000000;  % ppm atomic
Fsamp.Th230_Th232_2sig = 2 * Fsamp.Th230_Th232 * ...
                         sqrt((sampTh_blnk.Th232_Th229_err/sampTh_blnk.Th232_Th229)^2 + ...
                         (sampTh_blnk.Th230_Th229_err/sampTh_blnk.Th230_Th229)^2);

%% --------------------------------------------------------------------- %%
% Put it into Age Calc.
[Fages]          = MIT_Age_calc(Fsamp); 
dateOfChem       = datevec(powder.date_of_chem);
Fages.cor_BP     = Fages.cor - (dateOfChem(1)-1950);  % Year of chemistry - 1950
Fages.cor_BP_err = max(Fages.cor_2sig, Fages.uncor_2sig);

%% --------------------------------------------------------------------- %%
% Calculate impact of blank correction on final 238/236, 234/238, 230/229 
% ratios (e.g. (1-AD5/U5)*100%))
blank_impact.U238_U236   = (1 - sampU_blnk.U238_U236/sampU_spk.U238_U236)*100;
blank_impact.U238_U234   = (1 - sampU_blnk.U238_U234/sampU_spk.U238_U234)*100;
blank_impact.Th230_Th229 = (1 - sampTh_blnk.Th230_Th229/sampTh_spk.Th230_Th229)*100;
blank_impact.Th232_Th229 = (1 - sampTh_blnk.Th232_Th229/sampTh_spk.Th232_Th229)*100;
end