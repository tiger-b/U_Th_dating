%% ---------------------- MIT_tailyieldmassbias.m --------------------- %%
% DESCRIPTION
%   Applies tail, yield, and mass bias corrections to a sample or blank.
%
% INPUTS
%   stanTh1:  struct containing all raw Th data for first standard
%   stanTh2:  struct containing all raw Th data for second standard
%   sampTh:   struct containing all raw Th data for sample 
%   stanU1:   struct containing all raw U data for first standard
%   stanU2:   struct containing all raw U data for second standard
%   sampU:    struct containing all raw U data for sample 
%
% OUTPUTS
%   Fsamp:      struct of U-Th concentrations for reporting
%   Fages:      struct of final U-Th ages for reporting
%   blank_proc: struct of procedural blank information
%
% CONSTANTS
%   CONST:      struct of constants used for analysis (e.g., decay
%                  constants, spike U/Th ratios)
% 
% DEPENDENCIES
%   MIT_U_tail_corr_beta.m  :   Part 2 of U tail correction
% ----------------------------------------------------------------------- %
% AUTHORS: Christine Y Chen [CYC], David McGee [DM]
% Created:      06/09/15
% Last updated: 08/31/15
%% ---------------------- REVISION HISTORY ------------------------------ %
% 08/31/15    [CYC]     Changed name of 233.5 in tail correction to
%                       predicted and residual. 
%                       Changed U234 tail calc to not include 234.5 due to 
%                       issue of spuriously high 234.5 counts, matching the 
%                       08-17-2015 U-Th template. 
%                       Changed Th230/Th232 error and Th230/Th229 error to
%                       include U238 signal intensity, matching the 
%                       08-17-2015 U-Th template. 
%                       Added correction for U tailing in Th analysis based
%                       on measurement of 238U signal, matching the 
%                       08-17-2015 U-Th template. 
%                       Changed 234U tail estimate to be modified by residual 
%                       at 233.5 (even if it is negative), not maximum of 
%                       0.001 and 233.5 residual. 
% 08/13/15    [CYC]     Included Th 230 background into Th tail correction.
% 06/09/15    [CYC]     Created by CYC. Based on Excel reduction sheet from
%                       Dec 2014.
%
%% --------------------------------------------------------------------- %%
function [sampTh_yield, sampU_yield, U234t_234] = MIT_correctTailMassBiasYield(stanTh1, stanTh2, sampTh, ...
                                                                               stanU1,  stanU2,  sampU)

% Call global constants
global CONST;
%% ---------------------- Th tail correction --------------------------- %% 

stanTh1_tail = MIT_Th_tail_corr(stanTh1);
stanTh2_tail = MIT_Th_tail_corr(stanTh2);
sampTh_tail  = MIT_Th_tail_corr(sampTh);

    % ------------------------------------------------------------------- %
    % Sub-module that applies the tail correction to Th data.
    % ------------------------------------------------------------------- %
    function [cor] = MIT_Th_tail_corr(raw)
        % Load values of constants used in script
        global TAIL_Th230_Th232;
        global CPSPERVOLT;
        
        cor.Th230_Th232     = raw.Th230_Th232* ...
            (1 - TAIL_Th230_Th232*raw.Th232_V*CPSPERVOLT/raw.Th230_cps)*...
            (1 - raw.U238_V/raw.Th230_cps);       % Includes U tailing
        cor.Th230_Th232_err = sqrt(((cor.Th230_Th232 - raw.Th230_Th232)/4)^2 + ...
            raw.Th230_Th232_err^2 + ...
            ((raw.Th230_bg/raw.Th230_cps)* 0.25 * raw.Th230_Th232).^2 + ... % Includes Th230 background
            ((raw.U238_V/raw.Th230_cps)  * 0.10 * raw.Th230_Th232).^2);     % Includes raw U238 tailing
        cor.Th230_Th229     = raw.Th230_Th229* ...
            (1 - TAIL_Th230_Th232*raw.Th232_V*CPSPERVOLT/raw.Th230_cps)*...
            (1 - raw.U238_V/raw.Th230_cps);       % Includes U tailing;
        cor.Th230_Th229_err = sqrt(((cor.Th230_Th229 - raw.Th230_Th229)/4)^2 + ...
            raw.Th230_Th229_err^2 + ...
            ((raw.Th230_bg/raw.Th230_cps) * 0.25 * raw.Th230_Th229).^2 + ... % Includes Th230 background
            ((raw.U238_V/raw.Th230_cps)   * 0.10 * raw.Th230_Th229).^2);     % Includes raw U238 tailing
    end

%% ------------------ Th mass bias and yield correction ---------------- %%

% For Th230/Th232
sampTh_yield.Th230_Th232    = sampTh_tail.Th230_Th232 * ...
                              (CONST.std.MITh_1.Th230_Th229_mean/ ... 
                               CONST.std.MITh_1.Th232_Th229_mean)/...
                               (mean([stanTh1_tail.Th230_Th232, ...
                                      stanTh2_tail.Th230_Th232]));
sampTh_yield.Th230_Th232_err = sqrt((sampTh_tail.Th230_Th232_err^2 + ...
                                   ((stanTh1_tail.Th230_Th232 - stanTh2_tail.Th230_Th232)/...
                                     stanTh1_tail.Th230_Th232*sampTh_tail.Th230_Th232)^2));

% For Th232/Th229
sampTh_yield.Th232_Th229     = sampTh.Th232_Th229 * ...
                               CONST.std.MITh_1.Th232_Th229_mean/...
                               (mean([stanTh1.Th232_Th229, ...
                                      stanTh2.Th232_Th229]));
sampTh_yield.Th232_Th229_err = sqrt((sampTh.Th232_Th229_err^2 + ...
                                   ((stanTh1.Th232_Th229 - stanTh2.Th232_Th229)/...
                                     stanTh1.Th232_Th229*sampTh.Th232_Th229)^2));

% For Th230/Th229
sampTh_yield.Th230_Th229     = sampTh_tail.Th230_Th229 * ...
                               CONST.std.MITh_1.Th230_Th229_mean/...
                               (mean([stanTh1_tail.Th230_Th229, ...
                                      stanTh2_tail.Th230_Th229]));
sampTh_yield.Th230_Th229_err = sqrt((sampTh_tail.Th230_Th229_err^2 + ...
                                   ((stanTh1_tail.Th230_Th229-stanTh2_tail.Th230_Th229)/...
                                     stanTh1_tail.Th230_Th229*sampTh_tail.Th230_Th229)^2));

%% ---------------------- U tail correction ---------------------------- %%

% Send raw U ratios of standards and sample into two U tail correction
% sub-routines (MIT_U_tail_corr.m and MIT_U_tail_corr_beta.m).
stanU1_tail = MIT_U_tail_corr(stanU1); % Part 1
stanU2_tail = MIT_U_tail_corr(stanU2); % Part 1
sampU_tail  = MIT_U_tail_corr(sampU);  % Part 1

    % ------------------------------------------------------------------- %
    % Sub-module that applies the tail correction U data (Part 1).
    % ------------------------------------------------------------------- %
    function [cor] = MIT_U_tail_corr(raw)
        % Call the constants that will be used in this script 
        global U239_U238;        
        global CPSPERVOLT;  
        
        cor.U.U238_V        = raw.U.U238_V;
        cor.U.U236_V        = max([0.000001 raw.U.U236_V]);
        cor.U.U235_V        = cor.U.U238_V/raw.U.U238_U235;
        cor.U.U234_cps      = raw.U.U234_cps;
        cor.U.U233_V        = cor.U.U236_V/raw.U.U236_U233;
        
        cor.U.m             = (log(raw.Utail.U237_cps)-log(raw.Utail.U2365_cps))/...
                              (log(238-236.5));
        
        cor.Utail.U233_V    = exp(-cor.U.m*log(238-233) + log(raw.Utail.U237_cps)) + ...
            exp(-cor.U.m*log(236-233) + log(raw.Utail.U237_cps/cor.U.U238_V*cor.U.U236_V)) + ...
            exp(-cor.U.m*log(235-233) + log(raw.Utail.U237_cps/cor.U.U238_V*cor.U.U235_V));
        
        cor.Utail.U236_V    = exp(-cor.U.m*log(238-236) + log(raw.Utail.U237_cps));
        
        cor.Utail.U235_V    = exp(-cor.U.m*log(238-235) + log(raw.Utail.U237_cps)) + ...
            exp(-cor.U.m*log(236-235) + log(raw.Utail.U237_cps/cor.U.U238_V*cor.U.U236_V));
        
        cor.Utail.U2345_cps = exp(-cor.U.m*log(238-234.5) + log(raw.Utail.U237_cps)) + ...
            exp(-cor.U.m*log(236-234.5) + log(raw.Utail.U237_cps/cor.U.U238_V*cor.U.U236_V)) + ...
            exp(-cor.U.m*log(235-234.5) + log(raw.Utail.U237_cps/cor.U.U238_V*cor.U.U235_V));
        cor.U.U2345_cps     = raw.Utail.U2345_cps - cor.Utail.U2345_cps;
        
        cor.Utail.U2335_cps_predicted = exp(-cor.U.m*log(238-233.5) + log(raw.Utail.U237_cps)) + ...
            exp(-cor.U.m*log(236-233.5) + log(raw.Utail.U237_cps/cor.U.U238_V*cor.U.U236_V)) + ...
            exp(-cor.U.m*log(235-233.5) + log(raw.Utail.U237_cps/cor.U.U238_V*cor.U.U235_V));
        cor.Utail.U2335_cps_residual  = raw.Utail.U2335_cps - cor.Utail.U2335_cps_predicted;
        
        % Changed 234U tail estimate to be modified by residual at 233.5 
        % (even if it is negative), not maximum of 0.001 and 233.5 residual. 
        cor.Utail.U234_cps  = cor.U.U233_V*U239_U238*CPSPERVOLT + ...
            exp(-cor.U.m*log(238-234) + log(raw.Utail.U237_cps)) + ...
            exp(-cor.U.m*log(236-234) + log(raw.Utail.U237_cps/cor.U.U238_V*cor.U.U236_V)) + ...
            exp(-cor.U.m*log(235-234) + log(raw.Utail.U237_cps/cor.U.U238_V*cor.U.U235_V)) + ...
            cor.Utail.U2335_cps_residual;

% Old U234_cps (changed to above on 8/31/15, CYC)        
%         cor.Utail.U234_cps  = cor.U.U233_V*U239_U238*CPSPERVOLT + ...
%             exp(-cor.U.m*log(238-234) + log(raw.Utail.U237_cps)) + ...
%             exp(-cor.U.m*log(236-234) + log(raw.Utail.U237_cps/cor.U.U238_V*cor.U.U236_V)) + ...
%             exp(-cor.U.m*log(235-234) + log(raw.Utail.U237_cps/cor.U.U238_V*cor.U.U235_V)) + ...
%             exp(mean([log(max([0.001 cor.U.U2345_cps])) log(max([0.001 cor.U.U2335_cps]))]));
    end

% Part 2: "beta" correction
[stanU1_beta, sampU_beta, stanU2_beta] = MIT_U_tail_corr_beta(stanU1,      sampU,      stanU2, ...
                                                              stanU1_tail, sampU_tail, stanU2_tail);

%% ---------------------- U Yield correction --------------------------- %%  
% For U238/U234
sampU_yield.U238_U234     = sampU_beta.U.U238_U234_beta*CONST.std.CRM112a.U238_U234_mean/...
                            mean([stanU1_beta.U.U238_U234_beta, ...
                                  stanU2_beta.U.U238_U234_beta]);
sampU_yield.U238_U234_err = sampU_yield.U238_U234*...
                            sqrt((stanU1_beta.Utail.U238_U234_beta/stanU1_beta.U.U238_U234_beta/2)^2 + ...
                                 (stanU2_beta.Utail.U238_U234_beta/stanU2_beta.U.U238_U234_beta/2)^2 + ...
                                 (sampU_beta.Utail.U238_U234_beta/sampU_beta.U.U238_U234_beta)^2);
% For U238/U236
sampU_yield.U238_U236     = sampU_beta.U.U238_U236_beta;
sampU_yield.U238_U236_err = sampU_beta.Utail.U238_U236_beta;

%% ----------------------------- For final output ---------------------- %%
U234t_234 = sampU_beta.U.U234t_U234;
end