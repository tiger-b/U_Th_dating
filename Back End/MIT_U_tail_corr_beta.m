%% ------------------------ MIT_U_tail_corr_beta.m --------------------- %%
% DESCRIPTION
%   This module applies the tail correction U data (Part 2, "beta" 
%   correction).
%
% INPUTS
%   raw_stan1:  struct containing all raw U data for first standard
%   raw_samp:   struct containing all raw U data for sample 
%   raw_stan2:  struct containing all raw U data for second standard
%   tail_stan1: struct of U standard #1 data calculated during 1st tail corr
%   tail_samp:  struct of U sample data data calculated during 1st tail corr
%   tail_stan2: struct of U standard #2 data calculated during 1st tail corr
%
% OUTPUTS
%   beta_stan1: struct containing U data calculated during 2nd tail corr
%   beta_samp:  struct containing U data calculated during 2nd tail corr
%   beta_stan2: struct containing U data calculated during 2nd tail corr
%
% CONSTANTS
%   CPSPERVOLT, CONST
% ----------------------------------------------------------------------- %
% AUTHORS: Christine Y Chen [CYC], David McGee [DM]
% Created:      04/29/15
% Last updated: 08/14/15
%% ---------------------- REVISION HISTORY ------------------------------ %
% 
% 08/14/15    [CYC]     Changed 238/234 uncertainty related to mass bias 
%                       to be a function of 238/235 ratio for both samples
%                       and standards, according to DM's 6/15/15 changes 
%                       to the spreadsheet
% 05/01/15    [CYC]     Used global variable CONST to relay constants
%                       instead of having the function take the const 
%                       structure as an input.
% 04/29/15    [CYC]     Created by CYC. Based on Excel reduction sheet from
%                       Dec 2014.
%
%% --------------------------------------------------------------------- %%

function [beta_stan1, beta_samp, beta_stan2] = ...
         MIT_U_tail_corr_beta(raw_stan1,   raw_samp,  raw_stan2, ...
                              tail_stan1, tail_samp, tail_stan2)
% Call constants
global CPSPERVOLT;
global CONST;

% For standards, call only the beta correction for standards
beta_stan1 = U_tail_corr_beta_stan(raw_stan1, tail_stan1);
beta_stan2 = U_tail_corr_beta_stan(raw_stan2, tail_stan2);

% For samples, call both the beta correction for standards and samples.
beta_samp  = U_tail_corr_beta_stan(raw_samp, tail_samp);
beta_samp  = U_tail_corr_beta_samp(beta_samp, raw_samp, tail_samp, ...
                                   beta_stan1, beta_stan2);
    
    % ------------------------------------------------------------------- %
    % Sub-module that calculates beta corrected U values for standards
    % and samples. Standards need only call this function.
    % ------------------------------------------------------------------- %
    function [beta_stan] = U_tail_corr_beta_stan(raw_stan, tail_stan)
        beta_stan.U.U234t_U238     = tail_stan.Utail.U234_cps/tail_stan.U.U238_V/CPSPERVOLT;
        beta_stan.U.U234t_U234     = tail_stan.Utail.U234_cps/tail_stan.U.U234_cps;
        beta_stan.U.U235t_U235     = tail_stan.Utail.U235_V/tail_stan.U.U235_V/CPSPERVOLT;
        beta_stan.Utail.U234t_U234 = beta_stan.U.U234t_U234*0.1;
        beta_stan.Utail.U235t_U235 = beta_stan.U.U235t_U235*0.1;
        
        beta_stan.U.U238_U235      = raw_stan.U.U238_U235/...
                                     (1 - tail_stan.Utail.U235_V/CPSPERVOLT/tail_stan.U.U235_V);
        beta_stan.Utail.U238_U235  = sqrt((raw_stan.U.U238_U235_err/(1 - beta_stan.U.U235t_U235))^2 + ...
                                     (beta_stan.Utail.U235t_U235*raw_stan.U.U238_U235/(1 - beta_stan.U.U235t_U235)^2)^2);
        beta_stan.U.U238_U234      = raw_stan.U.U238_U234/(1 - beta_stan.U.U234t_U234);
        beta_stan.Utail.U238_U234  = sqrt((raw_stan.U.U238_U234_err/(1 - beta_stan.U.U234t_U234))^2 + ...
                                    (beta_stan.Utail.U234t_U234*raw_stan.U.U238_U234/(1 - beta_stan.U.U234t_U234)^2)^2);
        sqrt((raw_stan.U.U238_U234_err/(1 - beta_stan.U.U234t_U234))^2 + ...
                                     (beta_stan.Utail.U234t_U234*raw_stan.U.U238_U234/(1 - beta_stan.U.U234t_U234)^2)^2);
                                 
                                 
        beta_stan.U.beta           = (CONST.std.CRM112a.U238_U235_mean/beta_stan.U.U238_U235)^ ...
                                     (1/log(CONST.mass.U238/CONST.mass.U235));
        
        beta_stan.U.U238_U234_beta     = beta_stan.U.U238_U234*beta_stan.U.beta^...
                                         (log(CONST.mass.U238/CONST.mass.U234));
        beta_stan.Utail.U238_U234_beta = beta_stan.U.U238_U234_beta*...
                                         sqrt(4/3*(beta_stan.Utail.U238_U235/beta_stan.U.U238_U235)^2 + ...
                                         (beta_stan.Utail.U238_U234/beta_stan.U.U238_U234)^2);
    end

    % ------------------------------------------------------------------- %
    % Sub-module that calculates beta corrected U values for samples. 
    % Standards should call this function as well as the other function.
    % ------------------------------------------------------------------- % 
    function [beta_samp1] = U_tail_corr_beta_samp(beta_samp1, raw_samp, tail_samp, tail_stan_1, tail_stan_2)
        % For samples (must calculate for bracketing standards before calculating for sample)
        beta_samp1.U.U236t_U236     = tail_samp.Utail.U236_V/tail_samp.U.U236_V/CPSPERVOLT; 
        beta_samp1.Utail.U236t_U236 = beta_samp1.U.U236t_U236*0.1; 

        % Only calculated for samples
        beta_samp1.U.U238_U236      = raw_samp.U.U238_U236/...
                                      (1 - (tail_samp.Utail.U236_V/CPSPERVOLT)/tail_samp.U.U236_V); 
        beta_samp1.Utail.U238_U236  = sqrt((raw_samp.U.U238_U236_err/(1 - beta_samp1.U.U236t_U236))^2 + ...
                                      (beta_samp1.Utail.U236t_U236*raw_samp.U.U238_U236/(1 - beta_samp1.U.U236t_U236)^2)^2);

        % Only calculated for samples
        beta_samp1.U.U236_U233      = raw_samp.U.U236_U233*(1 - beta_samp1.U.U236t_U236);
        beta_samp1.Utail.U236_U233  = sqrt(raw_samp.U.U236_U233_err^2 + beta_samp1.Utail.U236t_U236^2);

        % Note that beta for the sample is calculated differently than beta
        % for the standards above. You take the average of two standard
        % values.
        beta_samp1.U.beta           = (CONST.std.CRM112a.U238_U235_mean/ ...
                                      mean([tail_stan_1.U.U238_U235 tail_stan_2.U.U238_U235]))^ ...
                                      (1/log(CONST.mass.U238/CONST.mass.U235));
                          
        % This calculation is also done by the standards, but with the beta
        % calculated for the sample.
        beta_samp1.U.U238_U234_beta     = beta_samp1.U.U238_U234*beta_samp1.U.beta^...
                                          (log(CONST.mass.U238/CONST.mass.U234));
        %beta_samp1.Utail.U238_U234_beta = beta_samp1.U.U238_U234_beta*...
        %                                  sqrt(4/3*(beta_samp1.Utail.U236_U233/beta_samp1.U.U236_U233)^2 + ...
        %                                  (beta_samp1.Utail.U238_U234/beta_samp1.U.U238_U234)^2); 
        beta_samp1.Utail.U238_U234_beta = beta_samp1.U.U238_U234_beta*...
                                          sqrt(4/3*(beta_samp1.Utail.U238_U235/beta_samp1.U.U238_U235)^2 + ...
                                          (beta_samp1.Utail.U238_U234/beta_samp1.U.U238_U234)^2);     
                                  
        % These are not calculated by the standards above
        beta_samp1.U.U238_U236_beta     = beta_samp1.U.U238_U236*beta_samp1.U.beta^...
                                          (log(CONST.mass.U238/CONST.mass.U236));
        beta_samp1.Utail.U238_U236_beta = beta_samp1.U.U238_U236_beta*...
                                          sqrt(2/3*(beta_samp1.Utail.U236_U233/beta_samp1.U.U236_U233)^2 + ...
                                          (beta_samp1.Utail.U238_U236/beta_samp1.U.U238_U236)^2);
        beta_samp1.U.U238_U235_beta     = beta_samp1.U.U238_U235*beta_samp1.U.beta^...
                                          (log(CONST.mass.U238/CONST.mass.U235));
        beta_samp1.Utail.U238_U235_beta = beta_samp1.U.U238_U235_beta*...
                                          sqrt((beta_samp1.Utail.U236_U233/beta_samp1.U.U236_U233)^2 + ...
                                          (beta_samp1.Utail.U238_U235/beta_samp1.U.U238_U235)^2);
        beta_samp1.U.U236_U233_beta     = beta_samp1.U.U236_U233*beta_samp1.U.beta^...
                                          (log(CONST.mass.U236/CONST.mass.U233));
        beta_samp1.Utail.U236_U233_beta = beta_samp1.U.U236_U233_beta*...
                                          sqrt((beta_samp1.Utail.U236_U233/beta_samp1.U.U236_U233)^2 + ...
                                          (beta_samp1.Utail.U236_U233/beta_samp1.U.U236_U233)^2);
                        
    end                   
end