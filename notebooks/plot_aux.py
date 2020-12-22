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
sys.path.append("/mnt/home/mrenzo/codes/python_stuff/plotFunc/")
from MESAreader import getSrcCol
import numpy as np

# ---------------------------------------------------------------------------------------
# Zeta Ophiuchi's data
# Spectral type: O9.5Vnn
# xi=15 # no idea what this is, maybe microturbulence


def get_zeta_oph_vsini():
    """
    returns the vsini and error from Zehe et al. 2019
    """
    vsini = 432  # km/s
    err_vsini = 16  # km/s
    return vsini, err_vsini


def zeta_oph_epsilon_he(ax):
    """
    plots the epsilon_he and error on ax
    data from Villamariz & Herrero 05
    """
    epsilon_he = 0.11  ## N(He/(N(H)+N(He)))
    err_epsilon_he_up = 0.05  # estimated from the range of values they explore
    err_epsilon_he_down = 0.02
    ax.errorbar(
        epslion_he, yerr=[err_epsilon_he_up, err_epsilon_he_down], fmt="o", color="r", zorder=1
    )


def zeta_oph_radius(ax):
    """
    plots the radius and error on ax
    data from Villamariz & Herrero 05
    """
    R = 8.3  # Rsun
    err_R = 1.5
    age = 5 # Myr -- this is based on fits from single rotating models, not trustworthy
    ax.errorbar(age, R, yerr=err_R, fmt="o", color="r", zorder=1)


def zeta_oph_mass(ax):
    """
    plots the mass and error on ax
    data from Villamariz & Herrero 05
    """
    M = 19  # Msun
    err_M = 11
    ax.errorbar(M, yerr=err_M, fmt="o", color="r", zorder=1)


def zeta_oph_spectroscopicHRD(ax):
    """
    plots Zeta ophiuchi logg and Teff on ax
    data from Villamariz & Herrero 05
    """
    logg = 3.7  # cm/s^2
    err_log_g = 0.15  # estimated from the range of values they explore
    Teff = 34000  # K
    err_Teff = 1500  # K # estimated from the range of Teff they explore
    # take log10
    log_Teff = np.log10(Teff)
    err_log_Teff = err_Teff / (Teff * np.log(10))
    ax.errorbar(log_Teff, logg, xerr=err_log_Teff, yerr=err_log_g, fmt="o", color="r", zorder=1)


def zeta_oph_HRD(ax):
    """
    plots Zeta ophiuchi L and Teff on ax
    data from Villamariz & Herrero 05
    """
    Teff = 34000  # K
    err_Teff = 1500  # K # estimated from the range of Teff they explore
    # take log10
    log_Teff = np.log10(Teff)
    err_log_Teff = err_Teff / (Teff * np.log(10))
    # visual magnitude
    # Mv = -4.2 # taken from Howarth & Prinja 1989
    L = 8.3e4  # Lsun
    err_L = 4e4
    # take log10
    log_L = np.log10(L)
    err_log_L = err_L / (L * np.log(10))
    ax.errorbar(log_Teff, log_L, xerr=err_log_Teff, yerr=err_log_L, fmt="o", color="r", zorder=10)
# -------------------------------------------------------------------
# radius plots
def get_radius_time(hfile):
    src, col = getSrcCol(hfile)
    t = src[:, col.index("star_age")] * 1e-6
    R = 10.0 ** (src[:, col.index("log_R")])
    return t, R


def plot_radius_time(ax, hfile1, c="#77CCCC", hfile2="", label=""):
    t, R = get_radius_time(hfile1)
    ax.plot(t, R, c=c, zorder=2, label=label)
    # post binary evolution part is optional
    if hfile2 != "":
        t, R = get_radius_time(hfile2)
        ax.plot(t, R, c=c, ls="-.", zorder=2, label=label)


# -------------------------------------------------------------------
# rotation
def get_surface_rotation_time(hfile):
    src, col = getSrcCol(hfile)
    t = src[:, col.index("star_age")] * 1e-6
    v = src[:, col.index("surf_avg_v_rot")]
    return t, v


def plot_surface_rotation_time(ax, hfile1, c="#77CCCC", hfile2="", label=""):
    t, v = get_surface_rotation_time(hfile1)
    ax.plot(t, v, c=c, zorder=2, label=label)
    # post binary evolution part is optional
    if hfile2 != "":
        t, v = get_surface_rotation_time(hfile2)
        ax.plot(t, v, c=c, zorder=2, ls="-.")


def get_omega_time(hfile):
    src, col = getSrcCol(hfile)
    omega_div_omega_crit = src[:, col.index("surf_avg_omega_div_omega_crit")]
    t = src[:, col.index("star_age")] * 1e-6  # in Myr
    return t, omega_div_omega_crit


def plot_omega_div_omega_crit(ax, hfile1, c="#77CCCC", ls="-", hfile2="", label=""):
    t, omega_div_omega_crit = get_omega_time(hfile1)
    ax.plot(t, omega_div_omega_crit, c=c, ls=ls, label=label, zorder=2)
    if hfile2 != "":
        t, omega_div_omega_crit = get_omega_time(hfile2)
        ax.plot(t, omega_div_omega_crit, c=c, ls="-.", zorder=2)


# -------------------------------------------------------------------
# Spectroscopic HRD
def getlogglogTeff(hfile):
    src, col = getSrcCol(hfile)
    log_Teff = src[:, col.index("log_Teff")]
    log_g = src[:, col.index("log_g")]
    return log_g, log_Teff


def plot_spectroscopic_HRD(ax, hfile1, c="#77CCCC", hfile2="", label=""):
    logg, logTeff = getlogglogTeff(hfile1)
    ax.plot(logTeff, logg, c=c, zorder=2, label=label)
    # post binary evolution part is optional
    if hfile2 != "":
        logg, logTeff = getlogglogTeff(hfile2)
        ax.plot(logTeff, logg, c=c, ls="-.", zorder=2)


# -------------------------------------------------------------------
# HRD
def getlogLlogTeff(hfile):
    src, col = getSrcCol(hfile)
    log_Teff = src[:, col.index("log_Teff")]
    logL = src[:, col.index("log_L")]
    return logL, log_Teff


def plot_HRD(ax, hfile1, c="#77CCCC", hfile2="", bfile="", label=""):
    logL, logTeff = getlogLlogTeff(hfile1)
    # if a binary history file is provided, print thicker lines during RLOF
    if bfile != "":
        src, col = getSrcCol(bfile)
        rl_relative_gap_1 = src[:, col.index("rl_relative_overflow_1")]
        ind = rl_relative_gap_1 > 0
        ax.plot(logTeff[ind], logL[ind], lw=8, c=c, zorder=2)
    # now make the plot
    ax.plot(logTeff, logL, c=c, zorder=2, label=label)
    # post binary evolution part is optional
    if hfile2 != "":
        logL, logTeff = getlogLlogTeff(hfile2)
        ax.plot(logTeff, logL, c=c, ls="-.", zorder=2)
# -------------------------------------------------------------------
# orbital velocity evolution
def plot_orbital_v2(ax, bfile, c="#77CCCC", label=""):
    srcb, colb = getSrcCol(bfile)
    v2 = srcb[:, colb.index("v_orb_2")]
    t = srcb[:, colb.index("age")]*1e-6
    ax.plot(t,v2, lw=3, c=c, label=label)
    ax.text(t[-1]+0.1, v2[-1], f"{v2[-1]:.1f}", fontsize=30)
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
    mp = 1.67262192369e-24  # grams
    mhe = 4 * mp  # neglects the binding energy of the he4 atom
    # get mass of the photosphere in grams
    mass = Mtot * Msun * dq
    # get number abundances from mass and mass fractions
    N_h = mass * X / mp
    N_he = mass * Y / mhe
    # calculate epsilon
    epsilon = N_he / (N_h + N - he)
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
    return 12 + np.log10(mass_frac / X)


def plot_surface_abundances(hfile1, hfile2="", ax="", label="", legend=False, do_log=True):
    """
    plot the surface abundances of a few isotopes
    the post binary evolution is optional
    """
    src, col = getSrcCol(hfile1)

    he4_1 = src[:, col.index("surface_he4")]
    h1_1 = src[:, col.index("surface_h1")]
    n14_1 = src[:, col.index("surface_n14")]
    c12_1 = src[:, col.index("surface_c12")]
    o16_1 = src[:, col.index("surface_o16")]
    t_1 = src[:, col.index("star_age")] * 1e-6

    if hfile2 != "":
        src2, col2 = getSrcCol(hfile2)
        he4_2 = src2[:, col2.index("surface_he4")]
        h1_2 = src2[:, col2.index("surface_h1")]
        n14_2 = src2[:, col2.index("surface_n14")]
        c12_2 = src2[:, col2.index("surface_c12")]
        o16_2 = src2[:, col2.index("surface_o16")]
        t_2 = src2[:, col2.index("star_age")] * 1e-6

        he4 = np.concatenate((he4_1, he4_2))
        h1 = np.concatenate((h1_1, h1_2))
        n14 = np.concatenate((n14_1, n14_2))
        c12 = np.concatenate((c12_1, c12_2))
        o16 = np.concatenate((o16_1, o16_2))
        t = np.concatenate((t_1, t_2))
    else:
        he4 = he4_1
        h1 = h1_1
        n14 = n14_1
        c12 = c12_1
        o16 = o16_1
        t = t_1

    if ax == "":
        # create plot figure if ax is not passed
        fig = plt.figure()
        gs = gridspec.GridSpec(100, 100)
        ax = fig.add_subplot(gs[:, :])

    if label != "":
        ax.set_title(label, fontsize=30)

    ax.plot(t, h1, c="b", label=r"$^1\mathrm{H}$")
    ax.plot(t, he4, c="r", label=r"$^4\mathrm{He}$")
    ax.plot(t, c12, c="g", label=r"$^{12}\mathrm{C}$")
    ax.plot(t, n14, c="m", label=r"$^{14}\mathrm{N}$")
    ax.plot(t, o16, c="y", label=r"$^{16}\mathrm{O}$")

    ax.axhline(h1[0], 0, 1, c="b", ls="--", lw=1)
    ax.axhline(he4[0], 0, 1, c="r", ls="--", lw=1)
    ax.axhline(n14[0], 0, 1, c="m", ls="--", lw=1)
    ax.axhline(c12[0], 0, 1, c="g", ls="--", lw=1)
    ax.axhline(o16[0], 0, 1, c="y", ls="--", lw=1)

    ax.set_xlabel(r"$\mathrm{time \ [Myr]}$")
    ax.set_ylabel(r"$\mathrm{Surface\ mass\ fraction}\ X_i$")
    # ax.set_xlim(xmin=8.5, xmax=10)
    if do_log:
        ax.set_yscale("log")
    if legend:
        ax.legend(ncol=2)
