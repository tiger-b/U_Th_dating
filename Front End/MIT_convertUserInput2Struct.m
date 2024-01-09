% ------------------------------------------------------------------- %
% Sub-module that creates a struct containing the powder weight, spike
% weight, spike type, solution volumes, and standard concentrations
% for a sample or a blank.
% ------------------------------------------------------------------- %
function [sampleInfo] = MIT_convertUserInput2Struct(pwderWt, spkWt, spkType, ...
    U_vol, Th_vol, ...
    std_U_conc, std_Th_conc, ...
    chemDate)
sampleInfo.samp_wt = pwderWt;
sampleInfo.spk_wt  = spkWt;
if strcmp(spkType, 'Calcite')
    sampleInfo.spk_type = 1;
else
    sampleInfo.spk_type = 0;
end
sampleInfo.U_vol        = U_vol;
sampleInfo.Th_vol       = Th_vol;
sampleInfo.std_U_conc   = std_U_conc;
sampleInfo.std_Th_conc  = std_Th_conc;
sampleInfo.date_of_chem = chemDate;
end
