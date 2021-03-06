# Instrument_processing
Codes for processing oceanographic instruments.

## Nortek_Signature

These are MATLAB function that I developed as part of a processing procedure for data collected from a set of Nortek Signature-500 ADCPs. In general, I've tried to provide "help" documentation/usage syntax within each function (though in some cases this may be incomplete). In my usage, these functions were employed together with a number of other project-specific functions (not included here) that do various corrections/quality control. A number of these codes are adapted from code that was provided by Nortek, and simplified or adjusted to suit my particular needs (particularly for beam mapping).

These codes were written based on variable names and data structures that are generated by converting the .ad2cp binaries to .mat format using Nortek's "MIDAS" software instead of the Signature Deployment software. In theory, the .mat files created with these two different methods contain all of the same data; however, variable name and structure format differ slightly between them. If you want to apply these codes to .mat files created with Signature Deployment software, it will be necessary to rename a number of variables but otherwise they should still work. I have the goal of eventually building a function that can "translate" between the two formats but haven't had a chance to yet.

Some of these functions are works-in-progress and may not produce accurate results (in particular, sigWavesProcess and sigBackscatter). 
