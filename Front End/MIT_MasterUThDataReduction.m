%% ------------------ MIT_MasterUThDataReduction.m --------------------- %%
% DESCRIPTION
%   Reads raw summary files exported from the mass spec, reduces data, and
%   outputs all relevant data into an Excel file. Prints messages to the 
%   command line indicating progress.
% ----------------------------------------------------------------------- %
% INPUTS
%   filename_sample_info: .csv file name of sample power information
%
% OUTPUTS
%   UTh-reduced-data-XX-XXX-XXX.xls : The relevant data for export. XX
%                                     represents the date of data
%                                     reduction.
%
% CONSTANTS
%   CONST:      struct of constants used for analysis (e.g., decay
%                  constants, spike U/Th ratios)
% 
% DEPENDENCIES
%   MIT_readRawSummaryData        : Reads raw summary data 
%   MIT_readConstants             : Loads pre-defined constants
%   MIT_convertBlank2Struct       : Converts blank data into struct format
%   MIT_printFlagMessages         : Prints flag messages for blank
%   MIT_convertUserInput2Struct   : Makes struct of blank powder info
%   MIT_calcFlagsAndYield         : Calculates percent chem yield and flags
%   MIT_correctTailMassBiasYield  : Tail, mass balance, and yield correction
%      MIT_Th_tail_corr           : Th tail correction
%      MIT_U_tail_corr            : Part 1 of U tail correction
%      MIT_U_tail_corr_beta       : Part 2 of U tail correction
%   MIT_convertSamp2Struct        : Converts sample data into struct format
%   MIT_reduceSample              : U-Th reduces one sample
%      MIT_spike_corr             : Spike correction
%      MIT_Age_calc               : Final age calculation
% ----------------------------------------------------------------------- %
% AUTHORS: Christine Y Chen [CYC], David McGee [DM], Ben Hardt [BH]
% Created:      04/30/15
% Last updated: 09/09/15
%
%% ---------------------- REVISION HISTORY ------------------------------ %
% 09/09/15   [CYC]     Renamed the function that takes user input and makes
%                      it into a struct 'MIT_convertUserInput2Struct'. It
%                      used to be 'MIT_makeSampInfoStruct'.
%                      Implemented BH's calculations for automation of
%                      significant figures. Note that upper and lower age
%                      significant figure calculations are not correct at
%                      the moment.
% 09/04/15   [CYC]     Edited a lot of stuff. Moved flags and yield
%                      calculation to its own function.
% 06/12/15   [CYC]     Created by CYC. Based on Excel reduction sheet from
%                      Dec 2014.
%
%% --------------------------------------------------------------------- %%
function MIT_MasterUThDataReduction(filename_sample_info)
versionNum = '1.1.090915';

fprintf('----------------- BEGIN U-TH DATA REDUCTION ----------------- \n')

% ------------------ Read in raw data from the csv files. --------------- %
fprintf('Reading raw data from summary files: ')
MIT_readRawSummaryData;
fprintf('[COMPLETE] \n')

% ---- Read in hard-coded constants (spikes, decay constants, masses) --- %
fprintf('Reading values of constants: ')
MIT_readConstants;
fprintf('[COMPLETE] \n')

% -- Declare global constants with set value throughout all scripts ----- %
global CPSPERVOLT;          % Counts per second in one volt
global U239_U238;           
global TAIL_Th230_Th232;    
global INIT230232;          % Initial 230Th/232Th ratio (in ppm atomic) 
global UNCERT230232;        % initial 230Th/232Th ratio 2-sigma uncertainty
global CONST;
CONST = constants;

% Set values to user-defined constant
TAIL_Th230_Th232 = tailTh230Th232;
U239_U238        = U239238; 
CPSPERVOLT       = countspersecpervolt;
INIT230232       = Th230Th232_ratio;
UNCERT230232     = Th230Th232_uncer;

% -------------- Convert blank values into a structure. ----------------- %
fprintf('Converting blank data into a structure: ')
MIT_convertBlank2Struct;
fprintf('[COMPLETE] \n')

% Put sample blank information into a structure.
blankInfo = MIT_convertUserInput2Struct(blank_pwderWt, blank_spkWt, blank_spkType, ...
                                        blank_U_vol, blank_Th_vol, ...
                                        blank_std_U_conc, blank_std_Th_conc, ...
                                        blank_dateOfChem);
                    
% ------------- Calculate percent yield and flags of blank -------------- %
fprintf('Calculating blank percent chemistry yield and testing flags:\n')
[blank_flags, blnk_yield] = MIT_calcFlagsAndYield(stdBlnkTh1_r, stdBlnkTh2_r, blankTh_r, ...
                                                  stdBlnkU1_r,  stdBlnkU2_r,  blankU_r, ...
                                                  blankInfo);
                                               
% Convert flag struct into an array and print any flag warning messages.                                               
blnk_flags = MIT_printFlagMessages(blank_flags);    
    
% ----- Apply tail, mass bias, and yield correction to the blank -------- %
fprintf('Tail, mass bias, and yield correction to the blank: ')
[blank.Th, blank.U, blank.U234t_234] = ...
    MIT_correctTailMassBiasYield(stdBlnkTh1_r, stdBlnkTh2_r, blankTh_r, ...
                                 stdBlnkU1_r,  stdBlnkU2_r,  blankU_r);
fprintf('[COMPLETE] \n')

%% ------------- Loop through each sample and reduce data -------------- %%
totalNumSamp = length(samp_name);
samp_flags = cell(numSamp, 6); % Pre-allocate space for flags
samp_yield =  NaN(numSamp, 2); % Pre-allocate space for percent chemistry yield
fprintf('\nReducing %d samples ... \n', totalNumSamp)
for i = 1:totalNumSamp
    
    % Find name of current sample being reduced
    sample_name = cellstr(samp_name(i));
    sample_name = sample_name{1, 1};
    fprintf(' \tReducing %s ... \n', sample_name)
    
    % -------- Read in powder weights, spike type and spike weight ------ %
    fprintf('\t\tReading in sample powder information:')
    sampInfo = MIT_convertUserInput2Struct(samp_pwderWt(i), samp_spkWt(i), samp_spkType(i), ...
                                           samp_U_vol(i), samp_Th_vol(i), ...
                                           samp_std_U_conc(i), samp_std_Th_conc(i), ...
                                           samp_dateOfChem(i));
    fprintf('[COMPLETE] \n')
    
    % ---------- Converts sample matrix data into a structure ----------- %
    fprintf('\t\tSample to structure conversion: ');
    [sampleTh_r, standardTh1_r, standardTh2_r, ...
     sampleU_r, standardU1_r, standardU2_r] = ...
                              MIT_convertSamp2Struct(i, ...
                                 Thsamp_data, Thstd1_data, Thstd2_data, ...
                                 Usamp_data, Usamp_tail_data, ...
                                 Ustd1_data, Ustd1_tail_data, ...
                                 Ustd2_data, Ustd2_tail_data);
    fprintf('[COMPLETE] \n')
    
    % ------------- Calculate percent yield and flags of blank -------------- %
    fprintf('\t\tCalculating blank percent chemistry yield and testing flags:\n')
    [sample_flags, samp_yield(i, :)] = MIT_calcFlagsAndYield(standardTh1_r, standardTh2_r, sampleTh_r, ...
                                                             standardU1_r,  standardU2_r,  sampleU_r, ...
                                                             sampInfo);
    
    % Convert flag struct into an array and print any flag warning messages.
    samp_flags(i, :) = MIT_printFlagMessages(sample_flags);
    
    % -------- Call data reduction script for this sample --------------- %
    fprintf('\t\tSample reduction: ')
    [Fconc(i), Fdates(i), Fblank_mass(i), Fblank_impact(i), U234t_U234(i)] = ...
        MIT_reduceSample(standardTh1_r, standardTh2_r, sampleTh_r, ...
        standardU1_r,  standardU2_r,  sampleU_r, ...
        sampInfo, blank);
    fprintf('[COMPLETE] \n')
end

fprintf('Reduction of %d samples: [COMPLETE] \n\n', totalNumSamp)

%% -------------- Generate output for the reduced data ----------------- %%
fprintf('Generating output for %d samples... \n', totalNumSamp)

% --------------- Collect some of data that you will output ------------- %

% Final concentrations
out_conc         = cell2mat(permute(struct2cell(Fconc),         [3 1 2])); 
out_dates        = cell2mat(permute(struct2cell(Fdates),        [3 1 2])); % U-Th Ages
out_blankMass    = cell2mat(permute(struct2cell(Fblank_mass),   [3 1 2])); % Blank masses
out_blankImpact  = cell2mat(permute(struct2cell(Fblank_impact), [3 1 2])); % Impact of blank on U-Th concentrations

% Re-arrange columns in out_dates for easier printing
out_dates = horzcat(out_dates(:, 1:6), out_dates(:, 9:10), out_dates(:, 7:8));

% 234tail/234 in per mil
out_samp_U234tU234 = U234t_U234;
out_blnk_U234tU234 = blank.U234t_234;

% Get sample 238, 236, 230, 229, 230background raw signal intensities
out_samp_U238raw    = Usamp_data (:, 9);
out_samp_U236raw    = Usamp_data (:, 11);
out_samp_Th230raw   = Thsamp_data(:, 9);
out_samp_Th229raw   = Thsamp_data(:, 11);
out_samp_Th230bgraw = Thsamp_data(:, 15);

% Get blank 238, 236, 230, 229, and 230background raw signal intensities
out_blnk_U238raw    = blankU_r.U.U238_V;
out_blnk_U236raw    = blankU_r.U.U236_V;
out_blnk_Th230raw   = blankTh_r.Th230_cps;
out_blnk_Th229raw   = blankTh_r.Th229_V;
out_blnk_Th230bgraw = blankTh_r.Th230_bg; 

% Get date of data reduction
t = datetime('now'); 
day = date;
out_dateTimeOfReduction = datestr(t); clear t;
out_dateOfReduction     = datestr(day); clear day;

%% --------------- Read everything to an external file ----------------- %%

out_filename = sprintf('UTh-reduced-data-%s.xls', out_dateOfReduction);
fprintf('Writing output to %s ... \n', out_filename)

% ----------- Write informational data to external file ----------------- %
fileID = fopen(out_filename, 'w');
fprintf(fileID, 'Date of reduction: %s\n', out_dateTimeOfReduction);
fprintf(fileID, 'Version of script: %s \n\n', versionNum);

% ------------------- Filenames of raw data files ----------------------- %
fprintf(fileID, 'FILENAMES OF RAW DATA FILES\n');
fprintf(fileID, 'User input of sample data\t%s\n', filename_sample_info); 
fprintf(fileID, 'U samp file\t%s\n',               filename_U_samp); 
fprintf(fileID, 'Th samp file\t%s\n',              filename_Th_samp); 
fprintf(fileID, 'U std file\t%s\n',                filename_U_std); 
fprintf(fileID, 'U tail file\t%s\n',               filename_U_tail); 
fprintf(fileID, 'Th std file\t%s\n\n',             filename_Th_std); 

% ------------- Export the values of various constants used ------------- %
fprintf(fileID, 'USER INPUT OF VALUES OF CONSTANTS\n');
fprintf(fileID, '230Th/232Th Tail Correction\t%0.4e\n', tailTh230Th232);
fprintf(fileID, 'U239/U238 ratio\t%0.4e\n', U239238);
fprintf(fileID, 'Detrital Th230/Th232 ratio\t%d\n', Th230Th232_ratio);
fprintf(fileID, 'Detrital Th230/Th232 uncertainty \t%0.4d\n\n', Th230Th232_uncer);

% ------- Export the sample data that the user originally inputted ------ %
fprintf(fileID, 'USER INPUT OF SAMPLE DATA\n');
fprintf(fileID, ['Sample Name\tSpike Type\tPowder Weight (g)\tSpike Weight (g)\t' ...
                 'U frac vol (mL)\tTh frac vol (mL)\t' ...
                 'Date of Chem\t' ...
                 'U std conc\tTh std conc\n']);

% Print out blank data first
blnk_name        = cellstr(blank_name);
blnk_name        = blnk_name{1, 1};
blnk_spk_type    = cellstr(blank_spkType);
blnk_spk_type    = blnk_spk_type{1, 1};
blnk_chem_date   = cellstr(blank_dateOfChem);
blnk_chem_date   = blnk_chem_date{1, 1};

% Print out blank name, spike type, and spike weight (skip powder weight)
fprintf(fileID, '%s\t%s\t\t%0.12f\t', ... 
        blnk_name, blnk_spk_type, blank_spkWt);
% Print out the blank volume of U and Th fractions
fprintf(fileID, '%0.2f\t%0.2f\t', ...
        blank_U_vol, blank_Th_vol);  
% Print out the blank date of chemistry    
fprintf(fileID, '%s\t', blnk_chem_date);
% Print out the blank U and Th standard concentrations used
fprintf(fileID, '%0.2f\t%0.2f\n', blank_std_U_conc, blank_std_Th_conc);
    
% Print out the sample data             
for i = 1:totalNumSamp
    % Find name of sample being written to output
    sample_name = cellstr(samp_name(i));
    sample_name = sample_name{1, 1};
    spk_type    =  cellstr(samp_spkType(1));
    spk_type    =  spk_type{1, 1};
    chem_date   = cellstr(samp_dateOfChem(i));
    chem_date   = chem_date{1, 1};
    
    % Print out sample name, powder weight, spike type and spike weight
    fprintf(fileID, '%s\t%s\t%0.12f\t%0.12f\t', ...
        sample_name, spk_type, samp_pwderWt(i), samp_spkWt(i));
    
    % Print out the sample volume of the U and Th fractions
    fprintf(fileID, '%0.2f\t%0.2f\t', ...
         samp_U_vol(i), samp_Th_vol(i));
    
    % Print out the sample date of chemistry
    fprintf(fileID, '%s\t', chem_date);
    
    % Print out the sample U and Th standard concentrations used
    fprintf(fileID, '%0.2f\t%0.2f\n', samp_std_U_conc(i), samp_std_Th_conc(i));
end

% ----------------------- Export analysis numbers ----------------------- %
fprintf(fileID, '\nANALYSIS NUMBERS AND DATE OF ANALYSIS\n');
fprintf(fileID, ['Sample Name\t' ...
                 'U samp #\tDate\tTh samp #\tDate\tU samp tail #\tDate\t' ...
                 'U std1 #\tDate\tU std1 tail #\tDate\t'...
                 'U std2 #\tDate\tU std2 tail #\tDate\t'...
                 'Th std1 #\tDate\tTh std2 #\tDate\n']);
for i = 1:totalNumSamp
    % Find name of sample being written to output
    if i == 1 % Print out name of procedural blank
        sample_name = blnk_name;
    else
        sample_name = cellstr(samp_name(i-1));
        sample_name = sample_name{1, 1};
    end
    
    % Find date and time of analysis to write to output
    samp_U_datetime      = cellstr(Usamp_datetime(i));
    samp_U_datetime      = samp_U_datetime{1, 1};
    samp_Th_datetime     = cellstr(Thsamp_datetime(i));
    samp_Th_datetime     = samp_Th_datetime{1, 1};
    samp_U_tail_datetime = cellstr(Usamp_tail_datetime(i));
    samp_U_tail_datetime = samp_U_tail_datetime{1, 1};
    std1_U_datetime      = cellstr(Ustd1_datetime(i));
    std1_U_datetime      = std1_U_datetime{1, 1};
    std1_U_tail_datetime = cellstr(Ustd1_tail_datetime(i));
    std1_U_tail_datetime = std1_U_tail_datetime{1, 1};
    std2_U_datetime      = cellstr(Ustd2_datetime(i));
    std2_U_datetime      = std2_U_datetime{1, 1};
    std2_U_tail_datetime = cellstr(Ustd2_tail_datetime(i));
    std2_U_tail_datetime = std2_U_tail_datetime{1, 1};
    std1_Th_datetime      = cellstr(Thstd1_datetime(i));
    std1_Th_datetime      = std1_Th_datetime{1, 1};
    std2_Th_datetime      = cellstr(Thstd2_datetime(i));
    std2_Th_datetime      = std2_Th_datetime{1, 1};

    fprintf(fileID, ['%s\t' ...
         '%d\t%s\t%d\t%s\t%d\t%s\t' ...
         '%d\t%s\t%d\t%s\t'...
         '%d\t%s\t%d\t%s\t'...
         '%d\t%s\t%d\t%s\n'], ...
         sample_name, ...
         samp_U_num(i), samp_U_datetime, samp_Th_num(i), samp_Th_datetime, samp_U_tail_num(i), samp_U_tail_datetime, ...
         std1_U_num(i), std1_U_datetime, std1_U_tail_num(i), std1_U_tail_datetime, ...
         std2_U_num(i), std2_U_datetime, std2_U_tail_num(i), std2_U_tail_datetime, ...
         std1_Th_num(i), std1_Th_datetime, std2_Th_num(i), std2_Th_datetime);
end

% ---- Export raw intensities, blank corr %, flags, and percent yield --- %
fprintf(fileID, '\nSAMPLE RAW INTENSITIES, IMPACT OF BLANK CORRECTION, STANDARD FLAGS, AND PERCENT YIELD\n');
fprintf(fileID, ['\t' ...
                 'Raw Intensity\tRaw Intensity\tRaw Intensity\tRaw Intensity\tRaw Intensity\tRaw Intensity\t' ...
                 'Blk Corr (%%)\tBlk Corr (%%)\tBlk Corr (%%)\tBlk Corr (%%)\t' ...
                 'Flag\tFlag\tFlag\tFlag\tFlag\tFlag\t' ...
                 'Yield\tYield\t\n']);
fprintf(fileID, ['Sample Name\t' ...
                 'U238 (V)\tU236 (V)\tTh230 (cps)\tTh229 (V)\tTh230bg (V)\t234tail/234 (per mil)\t' ...
                 '238/236\t238/234\t230/229\t232/229\t' ...
                 'Th mass bias\tU mass bias\tTh IC yield\tU IC yield\tTh low spike\tU low spike\t'...
                 'U236 (%%)\tTh229 (%%)\n']);

% First, output blank raw intensities, flags, and percent yield (skip impact of the blank correction)
fprintf(fileID, '%s\t', blnk_name);  % Name of blank
% Print out blank raw intensities (skip the impact of blank correction)
fprintf(fileID, '%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t\t\t\t\t', ...
         out_blnk_U238raw, out_blnk_U236raw, out_blnk_Th230raw, out_blnk_Th229raw, out_blnk_Th230bgraw, out_blnk_U234tU234);
% Print out blank flags
fprintf(fileID, '%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t', blnk_flags{:});
% Print out blank yield
fprintf(fileID, '%0.2f\t%0.2f\n', blnk_yield);

% Loop through and print the same information for the samples
for i = 1:totalNumSamp
    % Find name of sample being written to output
    sample_name = cellstr(samp_name(i));
    sample_name = sample_name{1, 1};
    
    % Print out sample name
    fprintf(fileID, '%s\t', ...
        sample_name);
    
    % Print out raw intensities
    fprintf(fileID, '%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t', ...
         out_samp_U238raw(i), out_samp_U236raw(i), out_samp_Th230raw(i), out_samp_Th229raw(i), out_samp_Th230bgraw(i), out_samp_U234tU234(i));
    
    % Print out impact of blank
    fprintf(fileID, '%0.2f\t%0.2f\t%0.2f\t%0.2f\t', out_blankImpact(i, :));

    % Print out flags
    fprintf(fileID, '%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t', samp_flags{i, :});
    
    % Print out percent yield for U236 and Th229
    fprintf(fileID, '%0.2f\t%0.2f\n', samp_yield(i, :));
end

% ------------------ Export blank mass data ----------------------------- %
fprintf(fileID, '\nPROCEDURAL BLANK MASSES\n');
fprintf(fileID, '230Th (fg)\t232Th (pg)\t234U (fg)\t238U (pg)\n');
fprintf(fileID, '%0.12f\t%0.12f\t%0.12f\t%0.12f\n', ...
                out_blankMass(1, :));

% -------- Export final sample concentration data and ages -------------- %
fprintf(fileID, '\nFINAL SAMPLE CONCENTRATIONS AND AGES\n');
fprintf(fileID, ['Sample Name\t' ...
                 '238U\t±(2sig)\t232Th\t±(2sig)\td234U\t±(2sig)\t(230Th/238U)\t±(2sig)\t230Th/232Th\t±(2sig)\t'...
                 'Age (uncorr)\t±(2sig)\tAge (corr)\t±(2sig)\td234init\t±(2sig)\tAge (corr)\t±(2sig)\t' ...
                 'Upper Age (corr)\tLower Age (corr)\n']);
fprintf(fileID, ['\t' ...
                 'ng/g\t\tpg/g\t\tper mille\t\tactivity\t\tppm atomic\t\t'...
                 'yr\t\tyr\t\tper mille\t\tyr BP\t\t' ...
                 'yr\tyr\n']);
for i = 1:totalNumSamp
    % Find name of sample being written to output
    sample_name = cellstr(samp_name(i));
    sample_name = sample_name{1, 1};

    fprintf(fileID, ['%s\t' ...
         '%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t' ...
         '%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t'...
         '%0.12f\t%0.12f\n'], ...
         sample_name, out_conc(i, :), out_dates(i, :));
end

% ---- Export final sample concentration data and ages WITH SIG FIGS ---- %
fprintf(fileID, '\nFINAL SAMPLE CONCENTRATIONS AND AGES WITH SIGNIFICANT FIGURES\n');
fprintf(fileID, ['Sample Name\t' ...
                 '238U\t±(2sig)\t232Th\t±(2sig)\td234U\t±(2sig)\t(230Th/238U)\t±(2sig)\t230Th/232Th\t±(2sig)\t'...
                 'Age (uncorr)\t±(2sig)\tAge (corr)\t±(2sig)\td234init\t±(2sig)\tAge (corr)\t±(2sig)\t' ...
                 'Upper Age (corr)\tLower Age (corr)\n']);
fprintf(fileID, ['\t' ...
                 'ng/g\t\tpg/g\t\tper mille\t\tactivity\t\tppm atomic\t\t'...
                 'yr\t\tyr\t\tper mille\t\tyr BP\t\t' ...
                 'yr\tyr\n']);
for i = 1:totalNumSamp
    % Find name of sample being written to output
    sample_name = cellstr(samp_name(i));
    sample_name = sample_name{1, 1};
    
    % Calculate significant figures for concentrations
    for j = 1:length(out_conc)
        if mod(j, 2) == 0  % Determine sig figs for 2 sigma uncertainty
            numSigFig = 1 - (1 + floor(log10(abs(out_conc(i, j)/2))));
            out_conc_sigfig(i, j) = round(out_conc(i, j), numSigFig);
        else  % Determine sig figs for concentrations and ratios
            numSigFig = 1 - (1 + floor(log10(abs(out_conc(i, j+1)/2))));
            out_conc_sigfig(i, j) = round(out_conc(i, j), numSigFig);
        end
    end
    
    % Calculate significant figures for ages
    for j = 1:length(out_dates)
        if mod(j, 2) == 0  % Determine sig figs for 2 sigma uncertainty
            numSigFig = 1 - (1 + floor(log10(abs(out_dates(i, j)/2))));
            out_dates_sigfig(i, j) = round(out_dates(i, j), numSigFig);
        else  % Determine sig figs for ages
            numSigFig = 1 - (1 + floor(log10(abs(out_dates(i, j+1)/2))));
            out_dates_sigfig(i, j) = round(out_dates(i, j), numSigFig);
        end
    end

    fprintf(fileID, ['%s\t' ...
         '%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t' ...
         '%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t%0.12f\t'...
         '%0.12f\t%0.12f\n'], ...
         sample_name, out_conc_sigfig(i, :), out_dates_sigfig(i, :));
end

fclose(fileID);

fprintf('----------------- END U-TH DATA REDUCTION ----------------- \n')
end