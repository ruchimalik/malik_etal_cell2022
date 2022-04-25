Photometry data were acquired using	**TDT RZ5P Base Processor  and Synapse software package**.  
https://www.tdt.com/system/fiber-photometry-system/

Analysis was performed using custom scriots written in MATLAB

*The analysis codes includes the following:*
1) Tom Davidson's photometry analysis suite (**tjd-shared-code and tjd-shared-codes-extras**).  
These folders have to be added to MATLAB path

2) Custom scripts
  * **Processing_photoM.m** --this scripts loads raw photometry data, demodulates the signals, normalizes the demodulated signal, and saves the output for futher analysis
  * **Plotting_photoM.m** -- this script loads and plots the output from "Processing_photoM.m". Options in the scripts allow for visualization of behavioraly relevant       timewindows
 
  
  
 
