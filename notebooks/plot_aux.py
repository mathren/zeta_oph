## author: Mathieu Renzo

## Author: Mathieu Renzo <mrenzo@flatironinstitute.org>
## Keywords: files

## Copyright (C) 2020 Mathieu Renzo

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see http://www.gnu.org/licenses/.

import sys
# the line below import stuff from git@github.com:mathren90/plotFunc.git
sys.path.append('/mnt/home/mrenzo/codes/python_stuff/plotFunc/')
from MESAreader import *
import matplotlib.pyplot as plt
from plotDefaults import *


# -------------------------------------------------------------------
# radius plots
def get_radius_time(hfile):
    src, col = getSrcCol(hfile)
    t = src[:, col.index("star_age")]*1e-6
    R = 10.0**(src[:, col.index("logR")])
    return t,R


def plot_radius_time(ax, hfile1,c1='#77CCCC', hfile2="", label=""):
    t, R = get_radius_time(hfile1)
    ax.plot(t, R, c=c1, zorder=2, label=label)
    # post binary evolution part is optional
    if hfile2 != "":
        t, R = get_radius_time(hfile2)
        ax.plot(t, R, c=c1, ls='-.', zorder=2, label=label)
# -------------------------------------------------------------------
# rotation
def get_surface_rotation_time(hfile):
    src, col = getSrcCol(hfile)
    t = src[:, col.index("star_age")]*1e-6
    v = src[:, col.index("surf_avg_v_rot")]
    return t,v

def plot_surface_rotation_time(ax, hfile1, c1='#77CCCC', hfile2="", label=""):
    t,v = get_surface_rotation_time(hfile1)
    ax.plot(t, v, c=c1, zorder=2, label=label)
    # post binary evolution part is optional
    if hfile2 != "":
        t,v = get_surface_rotation_time(hfile2)
        ax.plot(t, v, c=c1, zorder=2, ls='-.')

def get_omega_time(hfile):
    src, col = getSrcCol(hfile)
    omega_div_omega_crit = src[:, col.index("surf_avg_omega_div_omega_crit")]
    t = src[:, col.index("star_age")]*1e-6 # in Myr
    return t, omega_div_omega_crit

def plot_omega_div_omega_crit(ax, hfile1, c1="#77CCCC", hfile2="", label=""):
    t, omega_div_omega_crit = get_omega_time(hfile1)
    ax.plot(t, omega_div_omega_crit, c=c1, label=label, zorder=2)
    if hfile2 != "":
        t, omega_div_omega_crit = get_omega_time(hfile2)
        ax.plot(t, omega_div_omega_crit, c=c1, ls='-.',zorder=2)

# -------------------------------------------------------------------
# Spectroscopic HRD
def getlogglogTeff(hfile):
    src, col = getSrcCol(hfile)
    log_Teff = src[:, col.index("log_Teff")]
    log_g = src[:, col.index("log_g")]
    return log_g, log_Teff


def plot_spectroscopic_HRD(ax, hfile1, c1='#77CCCC', hfile2="", label=""):
    logg, logTeff = getlogglogTeff(hfile1)
    ax.plot(logTeff, logg, c=c1, zorder=2, label=label)
    # post binary evolution part is optional
    if hfile2 != "":
        logg, logTeff = getlogglogTeff(hfile2)
        ax.plot(logTeff, logg, c=c1, ls='-.', zorder=2)
# -------------------------------------------------------------------
# HRD
def getlogLlogTeff(hfile):
    src, col = getSrcCol(hfile)
    log_Teff = src[:, col.index("log_Teff")]
    logL = src[:, col.index("log_L")]
    return logL, log_Teff


def plot_HRD(ax, hfile1, c='#77CCCC', hfile2="", bfile="", label=""):
    logL, logTeff = getlogLlogTeff(hfile1)
    # if a binary history file is provided, print thicker lines during RLOF
    if bfile != "":
        src, col = getSrcCol(bfile)
        rl_relative_gap_1 = src[:, col.index('rl_relative_overflow_1')]
        ind = rl_relative_gap_1 > 0
        ax.plot(logTeff[ind], logL[ind], lw=8, c=c, zorder=2)
    # now make the plot
    ax.plot(logTeff, logL, c=c, zorder=2, label=label)
    # post binary evolution part is optional
    if hfile2 != "":
        logL, logTeff = getlogLlogTeff(hfile2)
        ax.plot(logTeff, logL, c=c, ls='-.', zorder=2)
# -------------------------------------------------------------------
# chemical composition
def get_epsilon_he(X, Y, Mtot=20, dq=1e-8):
    """
    convert mass fractions from MESA to the epslion spectroscopists like

    inputs:

    X, Y mass fractions of hydrogen and helium respectively,
    Mtot total mass
    dq fraction of the total mass that is the surface, 1e-8 is the MESA default

    output: epsilon_he = N(He)/(N(H)+N(He)) with N number abundance (?)

    TODO: verify this makes sense
    """
    # From $MESA_DIR/const/public/const_def.f90
    mp = 1.67262192369e-24 # grams
    mhe = 4*mp # neglects the binding energy of the he4 atom
    # get mass of the photosphere in grams
    mass = Mtot*Msun*dq
    # get number abundances from mass and mass fractions
    N_h = mass*X/mp
    N_he = mass*Y/mhe
    # calculate epsilon
    epsilon = N_he/(N_h+N-he)
    return epsilon

def get_epsilon(mass_frac, X):
    """
    convert mass fractions from MESA to the epslion spectroscopists like

    inputs:

    X mass fraction of hydrogen
    mass_frac mass fractions of the wanted element

    output: 12+log(mass_frac/X)

    TODO: verify this makes sense
    """
    return 12+np.log10(mass_frac/X)
# -------------------------------------------------------------------
