from MESAreader import getSrcCol, Rsun_cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

def get_pressure_scale_height(pfile):
    """ assuming hydrostatic equilibrium H = P/(rho g)"""
    src, col = getSrcCol(pfile)
    P = 10.0**(src[:, col.index("logP")])
    rho = 10.0**(src[:, col.index("logRho")])
    g = 10.00**(src[:, col.index("log_g")])
    return P/(rho*g) # cm


def plot_rho_Hp_radius(pfile, X_onset=0.42):
    fig = plt.figure()
    gs = gridspec.GridSpec(10, 10)
    ax = fig.add_subplot(gs[:,:])
    ax2 = ax.twinx()
    src, col = getSrcCol(pfile)
    log_rho = src[:, col.index("logRho")]
    log_r = np.log10(src[:, col.index("radius")]*Rsun_cm)
    H = get_pressure_scale_height(pfile)
    logH = np.log10(H)
    h1 = src[:, col.index("h1")]
    icore = np.argmin(np.absolute(h1-X_onset))
    ax.plot(log_r, log_rho, c='b')
    ax2.plot(log_r, logH, color='r')
    ax.axvline(log_r[icore], 0,1, ls='--', lw=2, c='m')
    ax.set_xlabel(r"$\log_{10}(r/\mathrm{[cm]})$")
    ax.set_ylabel(r"$\log_{10}(\rho/\mathrm{[g\ cm^{-3}]})$", color='b')
    ax2.set_ylabel(r"$\log_{10}(H_p/ \mathrm{[cm]})$", color="red")



def plot_rho_Hp_mass(pfile, X_onset=0.42):
    fig = plt.figure()
    gs = gridspec.GridSpec(10, 10)
    ax = fig.add_subplot(gs[:,:])
    ax2 = ax.twinx()
    src, col = getSrcCol(pfile)
    log_rho = src[:, col.index("logRho")]
    m = src[:, col.index("mass")]
    H = get_pressure_scale_height(pfile)
    logH = np.log10(H)
    h1 = src[:, col.index("h1")]
    icore = np.argmin(np.absolute(h1-X_onset))
    ax.plot(m, log_rho, c='b')
    ax2.plot(m, logH, color='r')
    ax.axvline(m[icore], 0,1, ls='--', lw=2, c='m')
    ax.set_xlabel(r"$\mathrm{Mass}\ [M_\odot]$")
    ax.set_ylabel(r"$\log_{10}(\rho/\mathrm{[g\ cm^{-3}]})$", color='b')
    ax2.set_ylabel(r"$\log_{10}(H_p/ \mathrm{[cm]})$", color="red")
