#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import glob
import astropy.io.fits as fits
import FITS_tools
from FITS_tools import hcongrid

from astropy.table import Table
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA


# In[329]:


test = fits.open('spec-0267-51608-0576.fits')
test2 = fits.open('spec-0268-51633-0055.fits')


# In[330]:


test.info()


# In[4]:


test[1].header


# In[5]:


Table(test[2].data).colnames


# In[ ]:


#spectro flux ugriz


# In[6]:


print(test[1].data['loglam'])
print(test2[1].data['loglam'])


# In[189]:


get_ipython().run_cell_magic('time', '', "\n#read in all relevant data into lists\n\nspec_list = glob.glob('*spec*')\n\nspectra = []\ncolor_array = []\nlog_lam = []\nz = []\n\nfor i,spec in enumerate(spec_list):\n    openspec_coadd = fits.getdata(spec, extname = 'COADD')\n    spectra.append(Table(openspec_coadd)['flux'][0:2870].data)\n    log_lam.append(Table(openspec_coadd)['loglam'][0:2870].data)\n\n    openspec_specobj = fits.getdata(spec, extname = 'SPECOBJ')\n    color_array.append(Table(openspec_specobj)['SPECTROFLUX'].data[0])\n    z.append(Table(openspec_specobj)['Z'].data[0])")


# In[249]:


# convert to linear space, then convert to rest frame by dividing by 1+z

lam = [[]]*len(log_lam)
lam_rf =[[]]*len(log_lam)

for i in range(len(log_lam)):
    lam[i] = 10**log_lam[i]
    lam_rf[i] = lam[i]/(1+z[i])


# In[266]:


#finding the wavelength range where all spectra overlap

lim_min_ll = lam_rf[0][0]
for i in range(len(lam_rf)):
    if (lam_rf[i][0] > lim_min_ll):
        lim_min_ll = lam_rf[i][0]    
print('lim_min_loglam = ', lim_min_ll)

lim_max_ll = lam_rf[0][2869]
for i in range(len(lam_rf)):
    if (lam_rf[i][2869] < lim_max_ll):
        lim_max_ll = lam_rf[i][2869]
print('lim_max_loglam = ', lim_max_ll)


# In[272]:


#creating array of overlapping wavelength range

log_lam_mask = (lam_rf[0] >= lim_min_ll) & (lam_rf[0] <= lim_max_ll)
adj_lam_rf = lam_rf[0][log_lam_mask]
print(adj_lam_rf)
print(len(adj_lam_rf))


# In[302]:


print(adj_lam_rf[0])
print(lam_rf[1][1100])


# In[382]:


get_ipython().run_cell_magic('time', '', '\n#cutting spectra along overlapping region\n\nspectra_cut = [[0]]*(len(spectra))\n\nfor i in range (len(lam_rf)):\n    for j in range(len(lam_rf[i])):\n        if (lam_rf[i][j] >= adj_lam_rf[0]) & (spectra_cut[i][0] == 0):\n            spectra_cut[i] = spectra[i][j:(j + (len(adj_lam_rf)))] ')


# In[317]:


get_ipython().run_cell_magic('time', '', '\n#running PCA on the entire sample\n\npca = PCA(n_components = 1000)\n \npc = pca.fit(spectra_cut)    \ncomponents = pca.components_\nvar = pca.explained_variance_ \nvar_ratio = pca.explained_variance_ratio_\nsing_vals = pca.singular_values_\nmean = pca.mean_\n\n#append mean to components for plotting\nmean_and_comp = np.vstack([mean,components])')


# In[453]:


#plot the mean and components

plt.figure(figsize=(20,25))
labels = ['mean', 'component 1', 'component 2', 'component 3', 'component 4']
for i in range (5):
    plt.subplot(5,1,i+1)
    plt.plot(adj_lam_rf, mean_and_comp[i], label = labels[i])
    plt.legend(loc = 'upper left', fontsize = 'x-large')


# In[323]:


norm_var_ratio = (var_ratio).cumsum()/sum(var_ratio)
norm_sing_vals = sing_vals.cumsum()/sum(sing_vals)


# In[325]:


#scree plot
plt.figure(figsize = (5,5))
plt.subplot(2,1,1)
plt.loglog(var_ratio)

plt.subplot(2,1,2)
plt.semilogx(norm_var_ratio**2)
plt.xlabel('Eigenvalue Number')
plt.ylabel('Cumulative Eigenvalues')
plt.ylim(0.85, 1.01)
plt.xlim(1,1100)


# In[328]:


#------------------------------------------------------------
# Plot the sequence of reconstructions
fig = plt.figure(figsize=(15, 10))
fig.subplots_adjust(hspace=0, top=0.95, bottom=0.1, left=0.12, right=0.93)

for i, n in enumerate([0, 4, 8, 20]):
    ax = fig.add_subplot(411 + i)
    #ax.plot(wavelengths, spectra[0], '-', c='gray')
    ax.plot(adj_lam_rf, mean + np.dot(norm_var_ratio[0:n], components[0:n]), '-k')

    if i < 3:
        ax.xaxis.set_major_formatter(plt.NullFormatter())

    ax.set_ylim(-2, 25)
    ax.set_ylabel('flux')

    if n == 0:
        text = "mean"
    elif n == 1:
        text = "mean + 1 component\n"
        text += r"$(\sigma^2_{tot} = %.2f)$" % norm_var_ratio[n]
    else:
        text = "mean + %i components\n" % n
        text += r"$(\sigma^2_{tot} = %.2f)$" % norm_var_ratio[n]

    ax.text(0.02, 0.93, text, ha='left', va='top', transform=ax.transAxes)

fig.axes[-1].set_xlabel(r'${\rm wavelength\ (\AA)}$')
plt.show()


# In[ ]:


# color magnitude diagram


# In[334]:


#extracting red and green flux from color array

color_r = []
color_g = []
for i in range (len(color_array)):
    color_r.append(color_array[i][2])
    color_g.append(color_array[i][1]) 


# In[358]:


#converting flux to apparent mag
r_mag = -2.5*np.log10(color_r) + 22.5
g_mag = -2.5*np.log10(color_g) + 22.5


# In[466]:


#color-magnitude diagram to separate red and blue galaxies (ignoring green valley)

plt.figure(figsize = (7,7))
plt.scatter(np.array(r_mag), (np.array(g_mag) - np.array(r_mag)), s = 4, alpha = .3, color = 'navy')
plt.plot([16,21], [.95,.95], lw = 1, color = 'k')
plt.ylabel('g-r', fontsize = 'x-large')
plt.xlabel('r-mag [apparent]', fontsize = 'x-large')
plt.title('Color-Magnitude Diagram', fontsize = 'x-large')
plt.xlim(20.5,16.5)
plt.ylim(0,1.25)


# In[364]:


# create condition to separate red gals from blue

red_cond = np.array(g_mag) - np.array(r_mag) > .95


# In[391]:


#cut spectra sample into red and blue gals

spectra_cut_red_gals = [[]]*sum(red_cond == True)
spectra_cut_blue_gals = [[]]*sum(red_cond == False)
k = 0
j = 0

for i in range (len(red_cond)):
    if red_cond[i] == True:
        spectra_cut_red_gals[j] = spectra_cut[i]
        j = j+1
    else:
        spectra_cut_blue_gals[k] = spectra_cut[i]
        k = k+1


# In[402]:


get_ipython().run_cell_magic('time', '', '\n#run PCA on Red gals\n\npca_red = PCA(n_components = 1000)\n \n# do transformation\npc_red = pca_red.fit(spectra_cut_red_gals)\n    \ncomponents_red = pc_red.components_\nvar_red = pc_red.explained_variance_ \nvar_ratio_red = pc_red.explained_variance_ratio_\nsing_vals_red = pc_red.singular_values_\nmean_red = pc_red.mean_\n\n#append mean to components for plotting\nmean_and_comp_red = np.vstack([mean_red,components_red])')


# In[460]:


#plot the mean and components

plt.figure(figsize=(20,25))
for i in range (5):
    plt.subplot(5,1,i+1)
    plt.plot(adj_lam_rf, mean_and_comp_red[i], label = labels[i])
    plt.xticks(ticks = np.arange(4250, 6751,250))
    plt.legend(loc = 'upper left', fontsize = 'x-large')


# In[430]:


#finding the wavelengths of the emission features in components 2-4

print(adj_lam_rf[(mean_and_comp_red[2].tolist()).index(max(mean_and_comp_red[2]))], ' -- OIII?')
print(adj_lam_rf[(mean_and_comp_red[3].tolist()).index(max(mean_and_comp_red[3]))], ' -- H-gamma?')
print(adj_lam_rf[(mean_and_comp_red[4].tolist()).index(max(mean_and_comp_red[4]))], ' -- H-beta?')

# https://astrobites.org/wp-content/uploads/2018/03/Berg_Fig4_excerpt.png    


# In[404]:


get_ipython().run_cell_magic('time', '', '\n#running PCA on blue gals\n\npca_blue = PCA(n_components = 1000)\n \npc_blue = pca_blue.fit(spectra_cut_blue_gals)\n    \ncomponents_blue = pc_blue.components_\nvar_blue = pc_blue.explained_variance_ \nvar_ratio_blue = pc_blue.explained_variance_ratio_\nsing_vals_blue = pc_blue.singular_values_\nmean_blue = pc_blue.mean_\n\n#append mean to components for plotting\nmean_and_comp_blue = np.vstack([mean_blue,components_blue])')


# In[461]:


#plot the mean and components

plt.figure(figsize=(20,25))
for i in range (5):
    plt.subplot(5,1,i+1)
    plt.plot(adj_lam_rf, mean_and_comp_blue[i], label = labels[i])
    plt.legend(loc = 'upper left', fontsize = 'x-large')


# In[ ]:




