
# coding: utf-8

# In[17]:

from glob import glob
import os
from collections import namedtuple
import pickle

# numerical tools
import numpy as np
import scipy.stats
import pandas as pd

# plotting tools
from matplotlib import pyplot as plt

# interactive notebook features
from tqdm import tqdm as tqdm
from ipywidgets import interact

# meg analysis
import mne


# ## Find all available subjects
# 
# Define where you store your `camcan` data in the variable `camcanroot`.

# In[6]:

# camcanroot = '/Volumes/Seagate Expansion Drive/camcan'
# camcanroot = os.path.join('D:', 'camcan')
camcanroot = os.path.join('/data', 'group', 'FANS', 'camcan-meg', 'camcan165', 'camcan165')
megdataroot = os.path.join(camcanroot, 'cc700', 'mri', 'pipeline','release004', 'BIDSsep', 'megraw')
subjects = glob(os.path.join(megdataroot, 'sub-*'))
ids = [os.path.split(subject)[-1][4:] for subject in subjects]

print(f'{len(subjects)} subjects found in {megdataroot}')


# In[9]:

# filter out no-files
subjects = [subject for subject in subjects if os.path.isfile(subject + '/meg/rest_raw.fif')]
print(len(subjects))


# ## Find the demographic information
# 
# Read the demographic information from the .tsv file provided.

# In[10]:

subject_details = pd.DataFrame.from_csv(os.path.join(camcanroot, 'cc700-scored/participant_data.csv'))


# ## Loop over MEG data
# 
# 

# In[ ]:

sub_params = namedtuple('sub_params',
                        ['pid', 'slopes', 'age', 'gender',
                         'intercepts', 'rsquared'])

all_parameters = []
psds = []

for subject in tqdm(range(174)):
    # resting state file
    restfile = os.path.join(subjects[subject], 'meg/rest_raw.fif')
    # raw data
    raw = mne.io.read_raw_fif(restfile, verbose='WARNING')
    # pick gradiometers
    picks = mne.pick_types(raw.info, meg='grad', eeg=False, stim=False, eog=False, exclude='bads')
    # do the PSD analysis
    psd, freqs = mne.time_frequency.psd_welch(
        raw, picks=picks, fmin=2, fmax=24, tmin=1, tmax=601, n_fft=2000, n_overlap=1000,
        verbose='WARNING', n_jobs=4
    )
    # Do the linear regression
    findices = (freqs < 7) | (freqs > 14)
    linfits = [scipy.stats.linregress(freqs[findices], np.log10(psd.T[findices, grad]))
               for grad in range(psd.shape[0])]

    psds.append(psd)
    all_parameters.append(
        sub_params(pid=ids[subject],
                   slopes=[l.slope for l in linfits],
                   intercepts=[l.intercept for l in linfits],
                   rsquared=[l.rvalue**2 for l in linfits],
                   age=subject_details.loc[ids[subject]].age,
                   gender=subject_details.loc[ids[subject]].gender_code)
    )


# In[ ]:

# save data to file
with open('./pickles/psds.pickle', 'wb+') as f:
    pickle.dump(psds, f)
with open('./pickles/all_parameters.pickle', 'wb+') as f:
    pickle.dump(all_parameters, f)


# ### Average power-spectrum

# In[ ]:

pltdata = np.log10(np.stack([p.mean(axis=0) for p in psds], axis=-1))

# plt.figure()
# plt.fill_between(freqs,
#                  np.mean(pltdata, axis=-1)-scipy.stats.sem(pltdata, axis=-1),
#                  np.mean(pltdata, axis=-1)+scipy.stats.sem(pltdata, axis=-1))
# plt.plot(freqs, np.mean(pltdata, axis=-1), color='orange')
# plt.show()


# ### Boxplot of regression slope and average $r^2$ across electrodes for each subject

# In[ ]:

# y = np.array([p.slopes for p in all_parameters])

# plt.figure()
# plt.boxplot(y.T)
# plt.show()

y = np.array([p.rsquared for p in all_parameters])

# plt.figure()
# plt.boxplot(y.T)
# plt.show()


# ### Scatter plot of age vs 1/f noise

# In[ ]:

x = [p.age for p in all_parameters]
y = [np.mean(p.slopes) for p in all_parameters]

# plt.figure()
# plt.scatter(x, y)
# plt.show()

