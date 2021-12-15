#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 14:53:08 2021

@author: rcruces
"""
# ----------------------------------------------------------------------------
#
# Database and group gradients
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
for k in ('100', '400', '1000'):
    atlas='schaefer-'+k
    print(atlas)

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
    annot_lh_fs5= nb.freesurfer.read_annot(micapipe + '/parcellations/lh.schaefer-' + k + '_mics.annot')
    Ndim = max(np.unique(annot_lh_fs5[0]))
    annot_rh_fs5= nb.freesurfer.read_annot(micapipe + '/parcellations/rh.schaefer-' + k + '_mics.annot')[0]+Ndim
    
    # replace with 0 the medial wall of the right labels
    annot_rh_fs5 = np.where(annot_rh_fs5==Ndim, 0, annot_rh_fs5) 
    
    # fsaverage5 labels
    labels_fs5 = np.concatenate((annot_lh_fs5[0], annot_rh_fs5), axis=0)
    
    # Read label for conte69
    labels_c69 = np.loadtxt(open(micapipe + '/parcellations/schaefer-' + k + '_conte69.csv'), dtype=np.int)
    
    # Mask of the medial wall on fsaverage 5
    mask_fs5 = labels_fs5 != 0
    
    # ----------------------------------------------------------------------------
    # LOAD ALL MATRICES
    # ----------------------------------------------------------------------------
    os.chdir("/Users/rcruces/tmp/micaConn/")
    
    # Within group mean gradient correlation
    DataSets=['MICS', 'CamCAN', 'MSC', 'audiopath', 'sudmex', 'EpiC', 'EpiC2']
    for dataset in DataSets:
        if os.path.exists('FC/' + dataset):
            # ----------------------------------------------------------------------------
            # FUNCTIONAL GRADIENTS
            # ----------------------------------------------------------------------------
            FC = fg.load_mtxs('FC/' + dataset + '/', Ndim, fg.load_fc, atlas)
            
            # Mean matrix across the z axis
            FCmean = np.mean(FC, axis=2)
            
            # Plot the log matrix
            corr_plot = plotting.plot_matrix(FCmean, figure=(10, 10), labels=None, cmap='Reds')
            corr_plot.figure.savefig(dataset +'-'+k+ '_FC_mtx.png')
            
            # Calculate the gradients
            FCgm = GradientMaps(n_components=10, random_state=None, approach='dm', kernel='normalized_angle')
            FCgm.fit(FCmean, sparsity=0.9)
            
            # Gradient plots
            fg.plot_grad(gmfit=FCgm, Color=['red'], screeTitle='Scree plot FC', s3DTitle='FC Gradients', s3Dlabels='', save=True, saveP=dataset + '_FC_')
            
            # # Map gradients to original parcels
            grad = [None] * 10
            for i, g in enumerate(FCgm.gradients_.T[0:10,:]):
                grad[i] = map_to_labels(g, labels_fs5,  fill=np.nan, mask=mask_fs5)
                plot_hemispheres(fs5_lh, fs5_rh, array_name=grad[i], cmap='RdYlBu_r', nan_color=(0, 0, 0, 1),
                                 zoom=1.3, size=(900, 750), embed_nb=True,
                                 color_bar='right', layout_style='grid', 
                                 label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                                 screenshot=True, filename=dataset +'-'+k+ '_FC-G'+str(i+1)+'.png')   
            
            # MPC group correlation
            FCrho = fg.group_corr(Mtx=FC, group_gm=FCgm, corrplots=False, gmplots=False, Title='FC ', S=0.9)
            
            # save to csv file
            np.savetxt(dataset +'-'+k+ '_FC-rho.csv', FCrho, delimiter=',')
            
        if os.path.exists('SC/' + dataset):
            # ----------------------------------------------------------------------------
            # STRUCTURAL GRADIENTS
            # ----------------------------------------------------------------------------
            SC = fg.load_mtxs('SC/' + dataset + '/', Ndim, fg.load_sc, atlas)
            
            # Mean matrix across the z axis
            SCmean = np.mean(SC, axis=2)
            
            # Plot the matrix
            corr_plot = plotting.plot_matrix(SCmean, figure=(10, 10), labels=None, cmap='Purples', vmin=0, vmax=10)
            corr_plot.figure.savefig(dataset +'-'+k+ '_SC_mtx.png')
            
            # SC Left hemi
            gm_SC_L = GradientMaps(n_components=10, random_state=None, approach='dm', kernel='normalized_angle')
            gm_SC_L.fit(SCmean[0:Ndim, 0:Ndim], sparsity=0)
            
            # SC Right hemi
            gm_SC_R = GradientMaps(n_components=10, alignment='procrustes', kernel='normalized_angle') # align right hemi to left hemi
            gm_SC_R.fit(SCmean[Ndim:Ndim*2, Ndim:Ndim*2], sparsity=0, reference=gm_SC_L.gradients_)
            
            # Gradient plots
            fg.plot_grad(gmfit=gm_SC_L, gmfit2=gm_SC_R, Color=['purple', 'slateblue'], screeTitle='Scree plot SC', s3DTitle='SC Gradients', s3Dlabels='', legendlabels=['Left SC', 'Right SC'], save=True, saveP=dataset + '_SC_')
            
            # Left and right gradients concatenated
            SC_gradients = np.concatenate((gm_SC_L.gradients_, gm_SC_R.aligned_), axis=0)
            
            # Map gradients to original parcels
            grad = [None] * 10
            for i, g in enumerate(SC_gradients.T[0:10,:]):
                grad[i] = map_to_labels(g, labels_fs5, fill=np.nan, mask=mask_fs5)
                plot_hemispheres(fs5_lh, fs5_rh, array_name=grad[i], cmap='RdYlBu_r', nan_color=(0, 0, 0, 1),
                                 zoom=1.3, size=(900, 750), embed_nb=True,
                                 color_bar='right', layout_style='grid', 
                                 label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                                 screenshot=True, filename=dataset +'-'+k+ '_SC-G'+str(i+1)+'.png')   
            
            # group correlation
            SCrho = fg.group_corr(Mtx=SC[0:Ndim, 0:Ndim, :], Mtx2=SC[Ndim:Ndim*2, Ndim:Ndim*2, :], 
                                group_gm=gm_SC_L, group_gm2=gm_SC_R, corrplots=False, gmplots=False, Title='SC ', S=0)
            
            # save to csv file
            np.savetxt(dataset +'-'+k+ '_SC-rho.csv', SCrho, delimiter=',')
    
            
        if os.path.exists('GD/' + dataset):
            # ----------------------------------------------------------------------------
            # GEODESIC DISTANCE GRADIENTS
            # ----------------------------------------------------------------------------
            GD = fg.load_mtxs('GD/' + dataset + '/', Ndim, fg.load_gd, atlas)
            
            # Mean matrix across the z axis
            GDmean = np.mean(GD, axis=2)
            
            # Plot the log matrix
            corr_plot = plotting.plot_matrix(GDmean[0:Ndim, 0:Ndim], figure=(10, 10), labels=None, cmap='Blues')
            corr_plot.figure.savefig(dataset +'-'+k+ '_GD_mtx.png')
            
            # GD Left hemi
            gm_GD_L = GradientMaps(n_components=10, random_state=None, approach='dm', kernel='normalized_angle')
            gm_GD_L.fit(GDmean[0:Ndim, 0:Ndim], sparsity=0.8)
            
            # GD Right hemi
            gm_GD_R = GradientMaps(n_components=10, alignment='procrustes', kernel='normalized_angle'); # align right hemi to left hemi
            gm_GD_R.fit(GDmean[Ndim:Ndim*2, Ndim:Ndim*2], sparsity=0.8, reference=gm_GD_L.gradients_)
             
            # Gradient plots
            fg.plot_grad(gmfit=gm_GD_L, gmfit2=gm_GD_R, Color=['dodgerblue', 'teal'], screeTitle='Scree plot GD', s3DTitle='GD Gradients', s3Dlabels='', legendlabels=['Left GD', 'Right GD'], save=True, saveP=dataset + '_GD_')
            
            # Left and right gradients concatenated
            GD_gradients = np.concatenate((gm_GD_L.gradients_, gm_GD_R.aligned_), axis=0)
            
            # Map gradients to original parcels
            grad = [None] * 10
            for i, g in enumerate(GD_gradients.T[0:10,:]):
                grad[i] = map_to_labels(g, labels_fs5, fill=np.nan, mask=mask_fs5)
                plot_hemispheres(fs5_lh, fs5_rh, array_name=grad[i], cmap='RdYlBu_r', nan_color=(0, 0, 0, 1),
                                      zoom=1.3, size=(900, 750), embed_nb=True,
                                      color_bar='right', layout_style='grid', 
                                      label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                                      screenshot=True, filename=dataset +'-'+k+ '_GD-G'+str(i+1)+'.png')
            
            # group correlation
            GDrho = fg.group_corr(Mtx=GD[0:Ndim, 0:Ndim, :], Mtx2=GD[Ndim:Ndim*2, Ndim:Ndim*2, :], 
                                group_gm=gm_GD_L, group_gm2=gm_GD_R, corrplots=False, gmplots=False, Title='GD', S=0.8)
            
            # save to csv file
            np.savetxt(dataset +'-'+k+ '_GD-rho.csv', GDrho, delimiter=',')
            
        if os.path.exists('MPC/' + dataset):
            # ----------------------------------------------------------------------------
            # MPC GRADIENTS
            # ----------------------------------------------------------------------------
            MPC = fg.load_mtxs('MPC/' + dataset + '/', Ndim, fg.load_mpc, atlas)
            # Mean matrix across the z axis
            MPCmean = np.mean(MPC, axis=2)
            
            # Plot the log matrix
            corr_plot = plotting.plot_matrix(MPCmean, figure=(10, 10), labels=None, cmap='Greens')
            corr_plot.figure.savefig(dataset +'-'+k+ '_MPC_mtx.png')
            
            # Calculate the gradients
            MPCgm = GradientMaps(n_components=10, random_state=None, approach='dm', kernel='normalized_angle')
            MPCgm.fit(MPCmean, sparsity=0.9)
            
            # Gradient plots
            fg.plot_grad(gmfit=MPCgm, Color=['green'], screeTitle='Scree plot MPC', s3DTitle='MPC Gradients', s3Dlabels='', save=True, saveP=dataset + '_MPC_')
            
            # Map gradients to original parcels
            grad = [None] * 10
            for i, g in enumerate(MPCgm.gradients_.T[0:10,:]):
                grad[i] = map_to_labels(g, labels_fs5, fill=np.nan, mask=mask_fs5)
                plot_hemispheres(fs5_lh, fs5_rh, array_name=grad[i], cmap='RdYlBu_r', nan_color=(0, 0, 0, 1),
                                      zoom=1.3, size=(900, 750), embed_nb=True,
                                      color_bar='right', layout_style='grid', 
                                      label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                                      screenshot=True, filename=dataset +'-'+k+ '_MPC-G'+str(i+1)+'.png')    
            
            # MPC group correlation
            MPCrho = fg.group_corr(Mtx=MPC, group_gm=MPCgm, corrplots=False, gmplots=False, Title='MPC', S=0.9)
            
            # save to csv file
            np.savetxt(dataset +'-'+k+ '_MPC-rho.csv', MPCrho, delimiter=',')

