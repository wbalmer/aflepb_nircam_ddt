from __future__ import division

# =============================================================================
# IMPORTS
# =============================================================================

import os
import pdb
import sys

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

from spaceKLIP import database, coron1pipeline, coron2pipeline, coron3pipeline, pyklippipeline, imagetools, analysistools

# Set the input and output directories and grab the input FITS files.
reduced = False
cleanalign = False
aligned = False

pad = False
coadd = False
crop = False

if aligned and coadd:
    input_dir = './spaceklip/coadded/'
    fitsfiles = sorted([input_dir + f for f in os.listdir(input_dir) if f.endswith('calints.fits')])
elif aligned:
    # input_dir = './spaceklip/aligned/'
    input_dir = './spaceklip/aligned_e3/'
    fitsfiles = sorted([input_dir + f for f in os.listdir(input_dir) if f.endswith('calints.fits')])
elif cleanalign:
    input_dir = './spaceklip/nanreplaced/'
    fitsfiles = sorted([input_dir + f for f in os.listdir(input_dir) if f.endswith('calints.fits')])
elif reduced:
    input_dir = './spaceklip/stage2/'
    fitsfiles = sorted([input_dir + f for f in os.listdir(input_dir) if f.endswith('.fits')])
else:
    input_dir = './uncal_e3/'
    # input_dir = './uncal_e3/'
    fitsfiles = sorted([input_dir + f for f in os.listdir(input_dir) if f.endswith('.fits')])

output_dir = './spaceklip/'

# Initialize the spaceKLIP database and read the input FITS files.
Database = database.Database(output_dir=output_dir)
Database.read_jwst_s012_data(datapaths=fitsfiles,
                             psflibpaths=None,
                             bgpaths=None)

select_obs = [
              'JWST_NIRCAM_NRCA2_F200W_MASKRND_MASK335R_SUB320A335R',
              'JWST_NIRCAM_NRCALONG_F356W_MASKRND_MASK335R_SUB320A335R',
              'JWST_NIRCAM_NRCALONG_F444W_MASKRND_MASK335R_SUB320A335R',
              ]

Database.obs = {k:Database.obs[k] for k in select_obs}

if not reduced:
    coron1pipeline.run_obs(database=Database,
                           steps={'saturation': {'n_pix_grow_sat': 1,
                                                 'grow_diagonal': True},
                                  'refpix': {'odd_even_columns': True,
                                             'odd_even_rows': True,
                                             'nlower': 4,
                                             'nupper': 4,
                                             'nleft': 4,
                                             'nright': 4,
                                             'nrow_off': 0,
                                             'ncol_off': 0},
                                  'dark_current': {'skip': True},
                                  'jump': {'rejection_threshold': 4.,
                                           'three_group_rejection_threshold': 4.,
                                           'four_group_rejection_threshold': 4.},
                                  'ramp_fit': {'save_calibrated_ramp': False}},
                           subdir='stage1')
    
if not reduced:
    coron2pipeline.run_obs(database=Database,
                           steps={'outlier_detection': {'skip': False}},
                           subdir='stage2')

ImageTools = imagetools.ImageTools(Database)

if not cleanalign:
    ImageTools.update_nircam_centers()
    
if not cleanalign:
    ImageTools.subtract_median(types=['SCI', 'SCI_TA', 'SCI_BG', 'REF', 'REF_TA', 'REF_BG'],
                                   subdir='medsub')

if not cleanalign:
    ImageTools.fix_bad_pixels(method='bpclean+timemed+dqmed+medfilt',
                              bpclean_kwargs={'sigclip': 3,
                                              'shift_x': [-1, 0, 1],
                                              'shift_y': [-1, 0, 1]},
                              custom_kwargs={},
                              timemed_kwargs={},
                              dqmed_kwargs={'shift_x': [-1, 0, 1],
                                            'shift_y': [-1, 0, 1]},
                              medfilt_kwargs={'size': 3},
                              subdir='bpcleaned')
    
if not cleanalign:
    ImageTools.replace_nans(cval=0.,
                            types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                            subdir='nanreplaced')

if not aligned:
    ImageTools.blur_frames()

if not aligned:
    if pad:
        ImageTools.pad_frames(
                              npix=[32, 33, 32, 33],
                              # npix=[1, 144, 73, 72], # shortwave
                              types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                              cval=0.
                             )
        
if not aligned:
    ImageTools.recenter_frames(
        spectral_type='F8V',
    )
    
if not aligned:
    ImageTools.align_frames(
        subdir='aligned_e3'
    )
    
pyklippipeline.run_obs(database=Database,
                       kwargs={'mode': ['ADI+RDI'],
                               'annuli': [4],
                               # 'movement': [0.5],
                               'subsections': [3],
                               'numbasis': [1, 2, 3, 4, 5, 10, 20, 50, 100],
                               'algo': 'klip'},
                       subdir='klipsub_e3')

Analysis = analysistools.AnalysisTools(Database)

companions = [[.3173, .0663, 1e-4]]
# companions = [[.3147, .0531, 1e-4]]

starfile = './AFLepA.vot'

mstar_err = 0.0

Analysis.raw_contrast(starfile,spectral_type='F8V',companions=companions, plot_xlim=(0,3), subdir='rawcon_e3')
Analysis.raw_contrast(starfile,spectral_type='F8V',companions=companions, plot_xlim=(0,3), subdir='rawcon')

inj_seps = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0]

Analysis.calibrate_contrast(
                            companions=companions,
                            injection_seps=inj_seps,
                            plot_xlim=(0,3),
                            subdir='calcon_e3'
                           )

Analysis.extract_companions(companions, 
                            starfile, 
                            mstar_err, 
                            klmode=100,
                            spectral_type='F8V', 
                            fitmethod='nested',
                            fitkernel='diag',
                            subdir='companions_e3'
                           )