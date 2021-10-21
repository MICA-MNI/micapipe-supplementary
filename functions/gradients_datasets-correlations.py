#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 19:23:49 2021

@author: rcruces
"""
# ----------------------------------------------------------------------------
#
# Database gradients correlations
#
# Created by RRC on September 2021 (the second year of the pademic)
# ----------------------------------------------------------------------------

# Set the environment
import os
import numpy as np
import nibabel as nb
from brainspace.plotting import plot_hemispheres
from brainspace.mesh.mesh_io import read_surface
from brainspace.datasets import load_conte69
from brainspace.gradient import GradientMaps
from brainspace.utils.parcellation import map_to_labels
import matplotlib.pyplot as plt
from nilearn import plotting # only to plot the matrices
import scipy.stats as ss

# Set the working directory to load the helper functions
os.chdir("/Users/rcruces/git_here/micapipe-supplementary/functions")
import func_gradients as fg

# Here we define the atlas 
atlas='schaefer-400'

# ----------------------------------------------------------------------------
# SURFACES
# ----------------------------------------------------------------------------
# Set the directory where the connectomes are
os.chdir("/Users/rcruces/tmp/micaConn/micapipe_tutorials")

# Variables
micapipe='/Users/rcruces/git_here/micapipe'

# Load conte69
c69_lh, c69_rh = load_conte69()

# Load fsaverage5
fs5_lh = read_surface('freesurfer/fsaverage5/surf/lh.pial', itype='fs')
fs5_rh = read_surface('freesurfer/fsaverage5/surf/rh.pial', itype='fs')

# Load annotation file in fsaverage5
annot_lh_fs5= nb.freesurfer.read_annot(micapipe + '/parcellations/lh.schaefer-400_mics.annot')
Ndim = max(np.unique(annot_lh_fs5[0]))
annot_rh_fs5= nb.freesurfer.read_annot(micapipe + '/parcellations/rh.schaefer-400_mics.annot')[0]+Ndim

# replace with 0 the medial wall of the right labels
annot_rh_fs5 = np.where(annot_rh_fs5==Ndim, 0, annot_rh_fs5) 

# fsaverage5 labels
labels_fs5 = np.concatenate((annot_lh_fs5[0], annot_rh_fs5), axis=0)

# Read label for conte69
labels_c69 = np.loadtxt(open(micapipe + '/parcellations/schaefer-400_conte69.csv'), dtype=np.int)

# Mask of the medial wall on fsaverage 5
mask_fs5 = labels_fs5 != 0

# ----------------------------------------------------------------------------
# LOAD ALL MATRICES BY DATASET
# ----------------------------------------------------------------------------
# Set the working directory to the 'out' directory
os.chdir("/Users/rcruces/tmp/micaConn")

# Datasets
DataSets=['MICS', 'CamCAN', 'EpiC', 'EpiC2', 'audiopath', 'sudmex', 'MSC']
N = len(DataSets)

# Create empty matrices for the datasets
FCds = np.empty([Ndim*2, Ndim*2, N], dtype=float)
SCds = np.empty([Ndim*2, Ndim*2, N], dtype=float)  
GDds = np.empty([Ndim*2, Ndim*2, N], dtype=float)  
MPCds = np.empty([Ndim*2, Ndim*2, N], dtype=float)  

for i, dataset in enumerate(DataSets):
    
    # Datasets Functional Connectomes
    if os.path.exists('FC/' + dataset):
        FCds[:,:,i] = np.mean(fg.load_mtxs('FC/' + dataset + '/', Ndim, fg.load_fc),  axis=2)
        
    # Datasets Structural Connectomes
    if os.path.exists('SC/' + dataset):
        SCds[:,:,i] = np.mean(fg.load_mtxs('SC/' + dataset + '/', Ndim, fg.load_sc),  axis=2)

    # Datasets Geodesic Distance Connectomes
    if os.path.exists('GD/' + dataset):
        GDds[:,:,i] = np.mean(fg.load_mtxs('GD/' + dataset + '/', Ndim, fg.load_gd),  axis=2)

    # Datasets MPC Connectomes
    if os.path.exists('MPC/' + dataset):
        MPCds[:,:,i] = np.mean(fg.load_mtxs('MPC/' + dataset + '/', Ndim, fg.load_mpc),  axis=2)


# Group level correlations
fg.dataset_corr(Means=SCds, Title='SC', DataSets=DataSets)
fg.dataset_corr(Means=FCds, Title='FC', DataSets=DataSets)
fg.dataset_corr(Means=GDds, Title='GD', DataSets=DataSets)
fg.dataset_corr(Means=MPCds, Title='MPC', DataSets=DataSets)