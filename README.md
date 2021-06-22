# Colocalization
Colocalization algorithm written in Matlab, used for Epigenetic readers in Expansion Microscopy

# MainColocalization.m
This is the main code to be used to perform colocalization. This code only works for .tif files, please convert your data to this format
prior to using it.

The main parameters are as follow:

locRoi ==> A radius in pixel used around detected particle
chi2   ==> Certainty threshold for initial detection, the bigger the number the more certain the algorithm need to be to accept a detection
           of particle and thus the less particle it detects. It is however more likely that these detected particle are real.
FWHM   ==> Expected FWHM of the particles to be detected (in Pixel)

nsim   ==> Number of random simulations to be ran to calculate the p-value

filename ==> Filename to store the output (end with .xlsx, making an excel file with the results)

The code will first request the user to point it toward a cell staining file, currently looking only for .tif file which have DAPI in the name
This can be modified here in Line 26:
[fileName,folder,~] = uigetfile('*DAPI*.tif','Pick a tif file to analyze');

The code will load the cell staining data and segment the cell which will be later use to perform the following calculation only on the inside
of the nucleus of the cell and ignore the rest of the data. If this step is not needed the code can also be adapted if one does not want/need
to use a nucleus staining. 

The code will then ask the user to point to a protein file. The file will be loaded and the particles will be localized using 2D Gaussian 
fitting.

Finally, the code will ask the user to point to a DNA marker .tif file. Once again the data will be loaded. Since the DNA markers are not that
localized, rather a intensity continuum with part that are darker and part that are brighter, we use a different approach to standard colocalization.
Indeed, we essentially get the intensity from the DNA marker at places where we detected protein molecules, then we do a random simulation to check
what are the odds that the protein is sitting at that position out of coincidence (random) rather than from colocalization. Indeed, since the intensity
distribution of the DNA marker is rather inhomogeneous we could expect that random data point would still have decent colocalization spot. The simulation
allows us to check this and estimate the p-value for our colocalization. A p-value <0.05, in this case would mean that the match we have
between the protein localised and the DNA marker is not random.

# mainPCC.m
The mainPCC.m performs a very similar analysis but instead of comparing intensity and localized spot, performs a correlation between the protein channel
and the DNA marker channel, a high correlation would indicate that most of the bright protein spot are placed in very close proximity to the bright spots
on the DNA marker channel.

In the paper, we used this analysis to compare with our colocalization approach and in the current form, our approach performs better than the pearson 
correlation paper.




