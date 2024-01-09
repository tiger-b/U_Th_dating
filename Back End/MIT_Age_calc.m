%% -------------------------- MIT_Age_calc.m ------------------------- %%
% DESCRIPTION
%   Calculates final U-Th ages from processed U-Th concentrations.
%
% INPUTS
%   conc:  struct containing the proceesed U-Th concentrations for sample
%
% OUTPUTS
%   ages:  struct containing final U-Th age information
%
% CONSTANTS
%   INIT230232:   initial 230Th/232Th ratio (in ppm atomic)
%   UNCERT230232: 2-sigma uncertainty
%   CONST:        struct of other constants
% ----------------------------------------------------------------------- %
% AUTHORS: Christine Y Chen [CYC], David McGee [DM], Benjamin Hardt [BFH]
% Created:      04/30/15
% Last updated: 09/16/15
%% ---------------------- REVISION HISTORY ------------------------------ %
% 09/16/15    [CYC]  Temporary editted the print statement to not be
%                    specific about which sample was infinite in order for
%                    Ben to test the script.  Will introduce the print
%                    statement at a later date.
% 09/04/15    [CYC]  Added to ouput struct the upper and lower bound of the
%                    uncorrected and corrected ages.
% 08/17/15    [CYC]  Changed from symbolic solver (solve) to numerical 
%                    solver (vpasolve) to improve speed. Added try-catch 
%                    statement so that for-loop may proceed after finding 
%                    infinite age.
% 06/15/15     [DM]  Changed line 69 to divide d234u by 2; expression
%                    within "sqrt" calculates the 1-sigma uncertainty, and
%                    the d234u is 2-sigma.
% 05/01/15    [CYC]  Used global variable CONST to relay constants
%                    instead of having the function take the const 
%                    structure as an input.
% 05/01/15    [CYC]  Specified which variable the solve functions should
%                    attempt to solve for (t) in order to speed up the
%                    calculation.
% 05/01/15    [DM]   Added 0.5 factor to calculation of d234init
%                    uncertainty. 
% 04/29/15    [CYC]  Created by CYC. Based on old AgeCalc.m file from
%                    04/22/15, in which BFH changed "lines 48 and 49" 
%                    (in AgeCalc.m) to remove multiplication by 2 for
%                    d234u (d234u is already 2 sigma).
%% --------------------------------------------------------------------- %%
function [ages] = MIT_Age_calc(conc)
% Load values for initial 230Th/232Th ratio (in ppm atomic) and its 
% 2-sigma uncertainty
global INIT230232;
global UNCERT230232;
global CONST;   % Struct of constants

% Assign the decay constants and atomic masses to variables used in code
lam230 = CONST.decay.Th230; % 230Th decay constant
lam232 = CONST.decay.Th232; % 232Th decay constant
lam234 = CONST.decay.U234;  % 234U decay constant
lam238 = CONST.decay.U238;  % 238U decay constant
m238   = CONST.mass.U238;   % 238U atomic mass
m232   = CONST.mass.Th232;  % 232Th atomic mass

% Assign concentration values to variables used in code
d234  = conc.d234U;
d234u = conc.d234U_2sig; %leaves out halflife uncertainty
ar08  = conc.Th230_U238;
ar08u = conc.Th230_U238_2sig;
r02   = conc.Th230_Th232 * 1e-6;
r02u  = conc.Th230_Th232_2sig * 1e-6;

ar28=ar08/r02*lam232/lam230;
ar02init=INIT230232*1e-6*lam230/lam232;
ar02unc=ar02init*UNCERT230232/INIT230232;

try
    syms t
    age=vpasolve(ar08-(1-exp(-lam230*t)+(d234/1000)*(lam230/(lam230-lam234))*(1-exp((lam234-lam230)*t))), t);
    agehi1=vpasolve((ar08+ar08u)-(1-exp(-lam230*t)+(d234/1000)*(lam230/(lam230-lam234))*(1-exp((lam234-lam230)*t))), t);
    agelo1=vpasolve((ar08-ar08u)-(1-exp(-lam230*t)+(d234/1000)*(lam230/(lam230-lam234))*(1-exp((lam234-lam230)*t))), t);
    agelo2=vpasolve(ar08-(1-exp(-lam230*t)+((d234+d234u)/1000)*(lam230/(lam230-lam234))*(1-exp((lam234-lam230)*t))), t);
    agehi2=vpasolve(ar08-(1-exp(-lam230*t)+((d234-d234u)/1000)*(lam230/(lam230-lam234))*(1-exp((lam234-lam230)*t))), t);
    ageunchi=sqrt((agehi1-age)^2+(agehi2-age)^2);
    ageunclo=sqrt((agelo1-age)^2+(agelo2-age)^2);
    ageunc=(ageunchi+ageunclo)/2;
    %correcting for initial 230Th
    agecorr=vpasolve((ar08-ar28*ar02init*exp(-lam230*t))-(1-exp(-lam230*t)+(d234/1000)*(lam230/(lam230-lam234))*(1-exp((lam234-lam230)*t))), t);
    agehi1=vpasolve((ar08-ar28*(ar02init-ar02unc)*exp(-lam230*t))-(1-exp(-lam230*t)+(d234/1000)*(lam230/(lam230-lam234))*(1-exp((lam234-lam230)*t))),t);
    agelo1=vpasolve((ar08-ar28*(ar02init+ar02unc)*exp(-lam230*t))-(1-exp(-lam230*t)+(d234/1000)*(lam230/(lam230-lam234))*(1-exp((lam234-lam230)*t))),t);
    agelo2=vpasolve((ar08-ar28*(ar02init)*exp(-lam230*t))-(1-exp(-lam230*t)+((d234+d234u)/1000)*(lam230/(lam230-lam234))*(1-exp((lam234-lam230)*t))),t);
    agehi2=vpasolve((ar08-ar28*(ar02init)*exp(-lam230*t))-(1-exp(-lam230*t)+((d234-d234u)/1000)*(lam230/(lam230-lam234))*(1-exp((lam234-lam230)*t))),t);
    ageunchi=sqrt((agehi1-agecorr)^2+(agehi2-agecorr)^2);
    ageunclo=sqrt((agelo1-agecorr)^2+(agelo2-agecorr)^2);
    agecorrunc=(ageunchi+ageunclo)/2;
    
    d234init=d234*exp(lam234*agecorr);
    d234initu=2*sqrt((d234u/2*exp(lam234*agecorr))^2+(d234*exp(lam234*(agecorr+0.5*agecorrunc))-d234init)^2);
    
    % Output the final ages as a structure
    ages.uncor         = double(age);
    ages.uncor_2sig    = double(ageunc);
    ages.cor           = double(agecorr);
    ages.cor_2sig      = double(agecorrunc);
    ages.d234init      = double(d234init);
    ages.d234init_2sig = double(d234initu);
    
    % Add to the ages struct the high and low of the uncorrected and
    % corrected ages
    ages.cor_hi        = double(max([agehi1 agehi2]));
    ages.cor_lo        = double(min([agelo1 agelo2]));
catch ME
    fprintf('Cannot solve sample. Possible infinite age.\n');
end
end