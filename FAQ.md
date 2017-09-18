# FAST Frequently Asked Questions

If FAST does not work or crashes, you are most likely doing something wrong. Please, read the documentation in the parameter file in detail! More information can be found below. If - after reading all of this - the problem still occurs, please email me! 

**I try to run v1.0, but get the following message: “% RESTORE: Incompatible saved routine written by IDL version newer than this one not restored”**  
This means that the executable does not work with your idl version. You can either update to a newer idl version, or you can remake the executable with the older idl version. In the download directory you can find a program called “mkexe.idl”. Copy this program to the idl directory. Now type on the command line:
```
idl mkexe.idl -arg fast
```
Next, replace the old fast executable with the new executable and it should work.

**If you fit spectra, does FAST take into account the resolution?**  
Fast does not broaden the models to match the resolution. There are a few ways of dealing with this. One option is to first adjust the model libraries to correct resolution. Another way is to fit the spectra in bins that are larger than both the spectral and model resolution. The example spectra file shows how to do this. Just make sure that your wavelength increments within one bin are smaller than the spectral and model resolution. 

**Why is FAST so slow?**  
There are a few options that slow down the program significantly. In particular, writing out the best fits and / or the final cube, and running Monte Carlo simulations. Also the grid size determines the speed of the fitting. For a first fit it is best to start without all these options, and keep the step sizes of the input parameters large.

**I get the following error message “Unable to allocate memory to make array: Please consult the supplier of the application”, what do I do wrong?**  
This means you do not have enough memory for the grid you are trying to fit. Increase the step size of the input parameters or use a computer with more memory.

**The output stellar masses and SFR do not make any sense?**  
Check whether the zeropoint in the parameter file is indeed the zeropoint of the fluxes in your input catalog.

**After a certain object, the output catalog behaves very strange?**  
If an object has no detection in a certain band, you have to put the flux to -99. Note that FAST does not recognize -NaN and this can disturb the rest of the fitting.

**What are the output files made by FAST?**  
The main output file is [OUTPUT_FILE].fout. If you do not define OUTPUT_FILE, the name will be the same as the input photometric and / or spectroscopic catalogs. Additionally FAST makes an output directory for each new grid (note this will not be done if you use the same grid for a new catalog with the same filters), with the name ([LIBRARY]_ [RESOLUTION]_[IMF]_XXXX_XXXX). If you are only interested in the output fitting parameters, do not bother looking at this directory. However, if you write out the best fits and / or χ2 grids for each object, you can find them here. Also, you can find a copy of the output (.fout) catalog and the input parameter file, both with the same name as the input photometric and / or spectroscopic catalogs. 

**Where can I find the full χ2 grid?**  
If you set SAVE_CHI_GRID to 1, the χ2 grid will be saved in the grid directory ([LIBRARY]_[RESOLUTION]_[IMF]_ XXXX_XXXX) as “chi_[OUTPUT_FILE]_[ID].save”. You can read in this file in idl using the following command:
```
IDL> restore, chi_[OUTPUT_FILE]_[ID].save, /verb
```
This file contains the χ2 grid “chi” (dim: z, tau, metal, age, Av), the number of degrees of freedom “n_degrees”, the full stellar mass “mass” and star formation rate “SFR” grid, the χ2 threshold values “chi_thr”, the best-fit solutions “best” and their confidence intervals “low” and “high”, the best scaling parameter “b_scale”, and all input parameters “key”.

**How do I know which grid directory belongs to which grid?**  
The code names of the grid directories are rather cryptic, and look like [LIBRARY]_[RESOLUTION]_[IMF]_XXXX_XXXX. You can check the copy of the parameter file in the directory, or check the 00README files. Every time a new grid is made, the code name and corresponding grid information are written to 00README. 

**Why are the FAST and HYPERZ masses not consistent?**  
Because you are not doing a fair comparison. If you assume the same library, initial mass function (and I do not mean applying some correction factor), star formation history, etc. you should find the same answer. Also note that HYPERZ has no template error function, so in order to do a fair comparison, you have to set this option to zero.

**Why are the best-fit solution for some objects not enclosed by the 68% confidence intervals if I input photometric redshifts?**  
This is caused by the fact that EAZY (or any other photometric redshift code) and FAST use different template sets and thus produce different PDFs of z. While developing FAST we have explored many different methods to derive confidence intervals when EAZY output is used. Unfortunately there is no perfect solution, but the one that we choose is close and easy to implement in FAST (without writing wrappers around EAZY). However, for some galaxies the best-fit solution does not fall within the confidence intervals. 

**Why are the Charlot & Bruzual (2007/2008) models not included and how can I make new libraries for FAST?**  
Because the CB07 files are not officially published, I do not distribute them. However, it is very easy to make a library for these models or any other star formation history (SFH) you want to explore. Just make sure the naming of the files is the same as the ma05 and bc03 models. We assume an exponentially declining SFH without gas recycling and no dust attenuation and use csp_galaxev from Bruzual & Charlot (2003) to apply the SFHs. Tau should be in log. The minimum and maximum log (tau/yr) and the stepsize depend on the grid you want to fit (see the parameter file). FAST reads in the binaries, so do not extract the models for specific ages.
