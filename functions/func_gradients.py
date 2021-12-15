#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 15:29:13 2021

@author: rcruces
"""

# ----------------------------------------------------------------------------
# GROUP GRADIENTS - HELPER FUNCTIONS
# ----------------------------------------------------------------------------
# Set the environment
import os
import glob
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

def load_fc(File, Ndim):
    """Loads and process a functional connectome"""
    
    # load the matrix
    mtx_fs = np.loadtxt(File, dtype=np.float, delimiter=' ')
    print(str(mtx_fs.shape))
    # slice the matrix remove subcortical nodes and cerebellum
    FC = mtx_fs[49:, 49:]
    
    # Remove the medial wall
    FC = np.delete(np.delete(FC, Ndim, axis=0), Ndim, axis=1)
    # Fishcer transform
    FCz = np.arctanh(FC)

    # replace inf with 0
    FCz[~np.isfinite(FCz)] = 0
    
    # Mirror the matrix
    FCz = np.triu(FCz,1)+FCz.T
    
    return FCz

def load_sc(File, Ndim):
    """Loads and process a functional connectome"""
    
    # load the matrix
    mtx_sc = np.loadtxt(File, dtype=np.float, delimiter=' ')
    
    # Mirror the matrix
    mtx_sc = np.log(np.triu(mtx_sc,1)+mtx_sc.T)
    mtx_sc[np.isneginf(mtx_sc)] = 0
    
    # slice the matrix remove subcortical nodes and cerebellum
    SC = mtx_sc[49:, 49:]
    SC = np.delete(np.delete(SC, Ndim, axis=0), Ndim, axis=1)
    
    # replace 0 values with almost 0
    SC[SC==0] = np.finfo(float).eps
    
    return SC

def load_gd(File, Ndim):
    """Loads and process a functional connectome"""
    
    # load the matrix
    mtx_gd = np.loadtxt(File, dtype=np.float, delimiter=' ')
    
    # Remove the Mediall Wall
    mtx_gd = np.delete(np.delete(mtx_gd, 0, axis=0), 0, axis=1)
    GD = np.delete(np.delete(mtx_gd, Ndim, axis=0), Ndim, axis=1)
    
    return GD

def load_mpc(File, Ndim):
    """Loads and process a functional connectome"""
    
    # load the matrix
    mtx_mpc = np.loadtxt(File, dtype=np.float, delimiter=' ')
    
    # Mirror the matrix
    MPC = np.triu(mtx_mpc,1)+mtx_mpc.T
    
    # Remove the medial wall
    MPC = np.delete(np.delete(MPC, 0, axis=0), 0, axis=1)
    MPC = np.delete(np.delete(MPC, Ndim, axis=0), Ndim, axis=1)
    
    return MPC

def load_mtxs(Path, Ndim, func, pattern):
    # List all the files
    # files = os.listdir(Path)
    files = glob.glob( Path + '*' + pattern + '_*' )
    
    # Load all the matrices
    M=np.empty([Ndim*2, Ndim*2, len(files)], dtype=float)
    for i, f in enumerate(files):
        print(f)
        M[:,:,i] = func(f, Ndim)
        
    return M

def plot_grad(gmfit, gmfit2=None, Color=['red', 'blue'], screeTitle='Scree plot', s3DTitle='Gradients', s3Dlabels='', 
              legendlabels=['gmfit1', 'gmfit2'], scree=True, scatter=True, aligned=False, save=False, saveP='file'):
    if scree == True:
        # Scree plot
        fig, ax = plt.subplots(1, 1, figsize=(7, 5))
        ax.plot(range(gmfit.lambdas_.size), gmfit.lambdas_, color = Color[0], linestyle='--', marker='o')
        if gmfit2 != None:
            ax.plot(range(gmfit2.lambdas_.size), gmfit2.lambdas_, color = Color[1], linestyle='dotted', marker='v')
            ax.legend(legendlabels)
        ax.set_xlabel('Component Nb')
        ax.set_ylabel('Eigenvalue')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.suptitle(screeTitle)
        if save == True:
            plt.savefig(saveP + 'scree.png', dpi=72)
            plt.close(fig)
        
    # Plot the gradients
    if aligned == False:
        g1=gmfit.gradients_[:, 0]
        g2=gmfit.gradients_[:, 1]
        g3=gmfit.gradients_[:, 2]
    else:
        g1=gmfit.aligned_[:, 0]
        g2=gmfit.aligned_[:, 1]
        g3=gmfit.aligned_[:, 2]
        
    if gmfit2 != None:
            g1R=gmfit2.aligned_[:, 0]
            g2R=gmfit2.aligned_[:, 1]
            g3R=gmfit2.aligned_[:, 2]
    
    if scatter == True:
        # Creating figure
        fig = plt.subplots(1, 1, figsize = (7, 5))
        ax = plt.axes(projection ="3d")
        
        # Creating plot
        ax.scatter3D(g1, g2, g3, color = Color[0])
        if gmfit2 != None:
            ax.scatter3D(g1R, g2R, g3R, color = Color[1], marker='v')
            ax.legend(legendlabels)
        plt.title(s3DTitle)
        ax.set_xlabel('G1')
        ax.set_ylabel('G2')
        ax.set_zlabel('G3')
        
        if save == True:
            plt.savefig(saveP + 'gradients.png', dpi=72)
            plt.close()
        else:
            plt.show()

def group_corr(Mtx, group_gm, Mtx2=None, group_gm2=None, corrplots=False, gmplots=False, Title='', S=0.9, Color=['gray']):
    # Empty array
    M=np.empty([Mtx.shape[2], group_gm.n_components], dtype=float)  

    for sub in range(Mtx.shape[2]):
        # Subject gradients
        gm_sub = GradientMaps(n_components=10, alignment='procrustes', kernel='normalized_angle'); # align to mean Gradient
        gm_sub.fit(Mtx[:,:,sub], sparsity=S, reference=group_gm.gradients_)
        
        if type(Mtx2) is np.ndarray:
            # Subject gradients
            gm_sub2 = GradientMaps(n_components=10, alignment='procrustes', kernel='normalized_angle'); # align to mean Gradient
            gm_sub2.fit(Mtx2[:,:,sub], sparsity=S, reference=group_gm2.aligned_)
            Color=['gray', 'black']
        
        if gmplots != False:
            # Gradient plots
            if not group_gm2:
                plot_grad(gmfit=gm_sub, Color=Color, screeTitle='sub-'+str(Title), s3DTitle=Title+'Gradients', s3Dlabels='', scree=False)
            else:
                plot_grad(gmfit=gm_sub, gmfit2=gm_sub2, Color=Color, screeTitle='sub-'+str(Title), s3DTitle=Title+'Gradients', s3Dlabels='', scree=False, aligned=True)
        
        if not group_gm2:
            # Variables to correlate
            x = gm_sub.aligned_
            y = group_gm.gradients_
        else:
            x = np.concatenate((group_gm.gradients_, group_gm2.aligned_), axis=0)
            y = np.concatenate((gm_sub.aligned_, gm_sub2.aligned_), axis=0)
        
        for i in range(group_gm.n_components):
            M[sub,i] = ss.spearmanr(x[:,i], y[:,i]).correlation 
        
        if corrplots != False:
            # Correlation plots
            fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, sharex=True, figsize=(12, 3))
            ax0.scatter(x[:, 0], y[:, 0], color = 'gray', marker='o')
            ax0.set_xlabel('Subject')
            ax0.set_ylabel('Group mean')
            ax0.spines['right'].set_visible(False)
            ax0.spines['top'].set_visible(False)
            ax0.set_title('G1 rho - ' + str(round(M[sub,0],2)))
            
            ax1.scatter(x[:, 1], y[:, 1], color = 'gray', marker='o')
            ax1.set_xlabel('Subject')
            ax1.set_ylabel('Group mean')
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)
            ax1.set_title('G2 rho - ' + str(round(M[sub,1],2)))
            
            ax2.scatter(x[:, 2], y[:, 2], color = 'gray', marker='o')
            ax2.set_xlabel('Subject')
            ax2.set_ylabel('Group mean')
            ax2.spines['right'].set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.set_title('G3 rho - ' + str(round(M[sub,2],2)))
            fig.suptitle( Title + ' sub' + str(sub))
            
    return M

def dataset_corr(Means, S=0.9, Title='', DataSets=[]):
    N = Means.shape[2]
    # Empty array of correlations by Gradient
    Gcorr = np.zeros([N, N, 10], dtype=float)
    
    for i in range(0, N-1):
        for j in range(i+1, N):
            #print(str(i) +', '+ str(j)) # print index
            
            if np.sum(Means[:,:,j]) != 0.0 and np.sum(Means[:,:,i]) != 0.0:    
                # Reference Database mean gradients
                gmRef = GradientMaps(n_components=10, alignment='procrustes', kernel='normalized_angle'); # align to mean Gradient
                gmRef.fit(Means[:,:,i], sparsity=S)
                
                # Moving Database mean gradients
                gmMov = GradientMaps(n_components=10, alignment='procrustes', kernel='normalized_angle'); # align to mean Gradient
                gmMov.fit(Means[:,:,j], sparsity=S, reference=gmRef.gradients_)
                
                # gradients components
                ref = gmRef.gradients_
                mov = gmMov.aligned_
                
                #  Align Database gradients
                for k in range(gmRef.n_components):    
                    Gcorr[j,i,k] = ss.spearmanr(ref[:,k], mov[:,k]).correlation 
            else:
                Gcorr[j,i,:] = 0

    for g in [0,1,2,9]:
        plotting.plot_matrix(Gcorr[:,:,g], figure=(10, 10), cmap='Reds', labels=DataSets, title=Title+' - G'+str(g+1), vmax=1, vmin=0 )
    
    # return Gcorr