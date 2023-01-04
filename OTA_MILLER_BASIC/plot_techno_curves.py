# -*- coding: utf-8 -*-

"""
Technological curves

65nm CMOS technology (TSMC)

Sylvain Favresse, 2022
From Marco Gonzalez's work on MATLAB (2021)

"""

import numpy as np
import matplotlib.pyplot as plt
import settings

def plot_techno_curves(techno_structs, techno_names):
    
    fig,ax = plt.subplots(2,2,figsize=(12,12))
    
    for i in range(len(techno_structs)):
        
        struct = techno_structs[i]
        
        ax[0][0].semilogx(struct['IN'], struct['GMID'],label=techno_names[i],linewidth=2)
        ax[0][1].plot(struct['GMID'], abs(struct['VDSAT']),label=techno_names[i],linewidth=2)
        ax[1][0].semilogy(struct['VGS'], struct['IN'],label=techno_names[i],linewidth=2)
        ax[1][1].plot(struct['VGS'], struct['GMID'],label=techno_names[i],linewidth=2)

    for i in range(2):
        for j in range(2):
            ax[i][j].legend(framealpha=1,fancybox=False,edgecolor='k')
            ax[i][j].grid(True,which='major')
            ax[i][j].grid(True,which='minor',alpha=0.3)
            
    
    ax[0][0].set_xlim(left=1e-10)
    ax[0][0].set_ylim((0,35))
    ax[0][1].set_ylim((0,0.8))
    ax[0][1].set_xlim((0,35))
    ax[1][0].set_xlim((0,1.2))
    ax[1][0].set_ylim(bottom=1e-10)
    ax[1][1].set_xlim((0,1.2))
    ax[1][1].set_ylim((0,35))
    
    ax[0][0].set_xlabel(r"$I_{Dn}$ [A]",fontsize=14)
    ax[0][0].set_ylabel(r"$g_m/I_D$ [1/V]",fontsize=14)
    ax[0][1].set_xlabel(r"$g_m/I_D$ [1/V]",fontsize=14)
    ax[0][1].set_ylabel(r"$V_{DS,sat}$ [V]",fontsize=14)
    ax[1][0].set_xlabel(r"$V_{GS}$ [V]",fontsize=14)
    ax[1][0].set_ylabel(r"$I_{Dn}$ [A]",fontsize=14)
    ax[1][1].set_xlabel(r"$V_{GS}$ [V]",fontsize=14)
    ax[1][1].set_ylabel(r"$g_m/I_D$ [1/V]",fontsize=14)
    
