# astro6410

This code take a query from the SDSS data release 16, cut on redshift (though not necessarily required), and runs Principal Component Analysis on the spectra. 

To get the data ready to be analyzed, I first converted all of the wavelengths to the rest frame by dividing by (1 + z). I then had to find the wavelength interval where all of the spectra overlap. After doing that, I cut each spectra on this interval. I ran PCA and plotted the mean as well as the first 4 components. The emission features are still a bit blurred, likely due to not considering the errors in redshift - though there could be other sources I haven't considered. 

To create a fuller picture of the data set and to better utilize PCA, I looked at the morphologies of the galaxies. I extracted the color array from the data set, then further extracted the red and green flux. I converted flux to magnitude by doing -2.5*log(flux) + 22.5 (where 22.5 was the zeropoint provided to me). I then created a color magnitude diagram and roughly separated out red and blue galaxies by creating a cut on g-mag - r_mag > 0.95. Those above this threshold I labeled as red and those below, blue. I did not consider the green valley.

I then reran PCA on each sample and plotted their means and first 4 components again. The results were as expected, the blue galaxies had strong H-alpha emission (signaling star formation) while the red galaxies did not. I also inspected some of the emission features of each sample and attempted to guess which elements were present. This analysis could be improved by considering errors in redshift and color flux as well as making more precise cuts on red, blue, and green galaxies. 
