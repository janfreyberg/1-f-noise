
# coding: utf-8

# In[ ]:

# path operations
from glob import glob
import os
from pathlib import Path

# data format and storage
from collections import namedtuple
import pickle

# numerical tools
import numpy as np
import scipy.stats
import pandas as pd

# plotting tools
from matplotlib import pyplot as plt
# %matplotlib notebook

# interactive notebook features
from tqdm import tqdm as tqdm
# from ipywidgets import interact

# meg analysis
import mne


# ## Find all available subjects
#
# Define where you store your `camcan` data in the variable `camcanroot`.

# In[ ]:

# camcanroot = Path('/Volumes') / 'Seagate Expansion Drive' /'camcan'
# camcanroot = Path('D:') / 'camcan'
camcanroot = Path('/data') / 'project' / 'CAMCAN' / \
    'camcan-meg' / 'camcan165'
# camcanroot = Path('/Users') / 'jan' / 'Documents' / 'eeg-data' / 'camcan'

megdataroot = camcanroot / 'cc700' / 'mri' / \
    'pipeline' / 'release004' / 'BIDSsep' / 'megraw'
subjects = list(megdataroot.glob('sub-*'))
ids = [os.path.split(subject)[-1][4:] for subject in subjects]

print(f'{len(subjects)} subjects found in {megdataroot}')

# filter out no-files
subjects = [subject for subject in subjects if (
    subject / 'meg' / 'task_raw.fif').is_file()]
print(f'{len(subjects)} subjects have resting-state recordings.')

# find the demographic info
subject_details = pd.DataFrame.from_csv(
    camcanroot / 'cc700-scored' / 'participant_data.csv')
print(f'Found subject information on {len([pid for pid in ids if pid in subject_details.index])} subjects.')


# ## Set up MEG analysis variables

# In[ ]:

veog = 'EOG062'
heog = 'EOG061'
ecg = 'ECG063'
recording = 'rest'
fmin, fmax = 2, 24


# ## Loop over MEG data
#
# First, make a data structure we can put the data in

# In[ ]:

sub_params = namedtuple('sub_params',
                        ['pid', 'slopes', 'age', 'gender',
                         'intercepts', 'rsquared'])


# For the pre-processing, we'll also use Maxwell-filtering to correct for
# extraneous influence. These are the necessary files:

# In[ ]:

ctfile = megdataroot / 'ct_sparse.fif'
ssscal = megdataroot / 'sss_cal.dat'


# Then, populate this data structure for each subject

# In[ ]:

all_parameters = []
psds = []

for subject in tqdm(range(len(subjects))):
    if (Path('.') / 'pickles' /
        (ids[subject] + '-' + recording + '.pickle')).is_file():
        continue
    # resting state file
    restfile = subjects[subject] / 'meg' / (recording + '_raw.fif')
    # raw data
    try:
        raw = mne.io.read_raw_fif(restfile, verbose='WARNING')
    except:
        continue
    # crop
    try:
        raw = raw.crop(tmin=20, tmax=500)
    except:
        continue
    # resample
    raw = raw.load_data()
    raw = raw.resample(256, n_jobs=10)
    # filter
    # filter the MEG data (exclude line noise)
    raw = raw.filter(0.5, 30, picks=mne.pick_types(raw.info, meg=True),
                     n_jobs=10)
    # filter the EOG data
    raw = raw.filter(0.5, 15, picks=mne.pick_types(
        raw.info, meg=False, eog=True),
        n_jobs=10)
    # filter the ECG data
    raw = raw.filter(0.5, 15, picks=mne.pick_types(
        raw.info, meg=False, ecg=True),
        n_jobs=10)
    # maxwell-correction
    raw = mne.preprocessing.maxwell_filter(raw, cross_talk=str(ctfile), calibration=str(ssscal),
                                           st_duration=10, st_correlation=0.98)
    # pick gradiometers
    picks = mne.pick_types(raw.info, meg='grad', eeg=False,
                           stim=False, eog=False, exclude='bads')
    # run an ICA
    try:
        ica = mne.preprocessing.run_ica(raw, n_components=0.95,
                                        picks=picks,
                                        eog_ch=veog, ecg_ch=ecg)
        raw = ica.apply(raw, exclude=ica.exclude)
    except:
        continue

    # do the PSD analysis
    psd, freqs = mne.time_frequency.psd_welch(
        raw, picks=picks, fmin=fmin, fmax=fmax, n_fft=2000, n_overlap=1000,
        verbose='WARNING', n_jobs=10
    )

    # Do the linear regression
    findices = (freqs < 7) | (freqs > 14)
    linfits = [scipy.stats.linregress(freqs[findices], np.log10(psd.T[findices, grad]))
               for grad in range(psd.shape[0])]

    psds.append(psd)
    all_parameters.append(
        sub_params(pid=ids[subject],
                   slopes=np.array([l.slope for l in linfits]),
                   intercepts=np.array([l.intercept for l in linfits]),
                   rsquared=np.array([l.rvalue**2 for l in linfits]),
                   age=subject_details.loc[ids[subject]].age,
                   gender=subject_details.loc[ids[subject]].gender_code)
    )
    with open(Path('.') / 'pickles' / (ids[subject] + '-' + recording + '.pickle'), 'wb+') as f:
        pickle.dump((psd, sub_params), f)


# ### Save data to pickles

# In[ ]:

# save data to file
with open('./pickles/psds.pickle', 'wb+') as f:
    pickle.dump(psds, f)
with open('./pickles/all_parameters.pickle', 'wb+') as f:
    pickle.dump(all_parameters, f)
