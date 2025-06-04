# U_Th_dating
Contains code to calculate U-Th ages with data from the McGee Lab MC-ICP-MS. 
Originally written by Christine Chen, Christopher Kinsley, Ben Hardt, and David McGee


Instructions:
1. Create a directory containing the "Back End" and "Front End" folders.
2. In the same directory, fill out an input sheet with sample information (sample numbers, date of chemistry, etc.).
	-If a run deviated from the normal sequence in any way, sample, standard, and tail numbers must all be input manually. Otherwise, 	only U and Th sample numbers are required input and the software will find the standard and tail numbers automatically. 
3. In MATLAB, navigate to this directory and add it to your search path. 
4. Run MIT_MasterUThDataReduction('sample_input.csv') to generate an output file in the same directory.

Notes:
1. The 230Th/232Th Tail Correction needs to be calculated separately.
2. The sample weight for the procedural blank needs to be input as 0.


Here, a "sample_data" folder with summary files and the excel reduction sheet are provided with "sample_input.csv" if you want to quickly do a test run on your machine.