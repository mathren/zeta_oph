import sys
# the line below import stuff from git@github.com:mathren90/plotFunc.git
sys.path.append('/mnt/home/mrenzo/codes/python_stuff/plotFunc/')
from MESAreader import *
import matplotlib.pyplot as plt
from plotDefaults import *
# plot_aux defines many auxiliary functions to make plots
from plot_aux import *
set_plotDefaults()

def plot_Dmix(pfile, ax="", legend=False):
    if ax == "":
        fig = plt.figure(figsize=(10,10))
        gs = gridspec.GridSpec(100, 100)
        ax = fig.add_subplot(gs[:, :])
    age = get_age_from_profile(pfile)
    timestamp = f"${age:.5f}"+r"\,\mathrm{Myr}$"

    src, col = getSrcCol(pfile)
    mass = src[:, col.index("mass")]
    D_mix = 10**(src[:, col.index('log_D_mix')])
    D_mix_non_rotation = 10.0**(src[:, col.index('log_D_mix_non_rotation')])
    D_conv = 10.0**(src[:, col.index('log_D_conv')])
    D_semi = 10.0**(src[:, col.index('log_D_semi')])
    D_thrm = 10.0**(src[:, col.index('log_D_thrm')])

    D_ov = D_mix_non_rotation - D_thrm - D_semi
    D_rot = D_mix - D_mix_non_rotation

    # ax.scatter(mass, D_mix, lw=10, label="$\mathrm{Total}$", c="c")
    ax.scatter(mass, D_rot, lw=0, s=5, label="$\mathrm{Rotation}$", c='#77AADD')
    ax.fill_between(mass, D_rot, alpha=0.75, fc='#77AADD')

    ax.scatter(mass, D_ov, lw=0, s=5, label="$\mathrm{Overshooting}$", c='b')
    ax.fill_between(mass, D_ov, alpha=0.75, fc='b')

    ax.scatter(mass, D_conv, lw=0, s=5, label="$\mathrm{Convection}$", c="r")
    ax.fill_between(mass, D_conv, alpha=0.75, fc='r')

    # ax.scatter(mass, D_semi, lw=0, s=5, label="$\mathrm{Semiconvection}$", c="#DDDD77")
    # ax.fill_between(mass, D_semi, alpha=0.75, fc='#DDDD77')

    ax.scatter(mass, D_thrm, lw=0, s=5, label="$\mathrm{Thermohaline}$", c="#CC99BB")
    ax.fill_between(mass, D_thrm, alpha=0.75, fc='#CC99BB')

    ax.text(0.05,0.1, timestamp, fontsize=30, va="center",
            transform=ax.transAxes, bbox=dict(facecolor='w',
                                              edgecolor='black', boxstyle='round,pad=0.2', alpha=0.9))
    ax.set_yscale('log')
    if legend:
        conv, = ax.plot(np.nan, np.nan, lw=30, ls='-', c='r',label="$\mathrm{Convection}$")
        ov, = ax.plot(np.nan, np.nan, lw=30, ls='-', c='b',label="$\mathrm{Overshooting}$")
        rot, = ax.plot(np.nan, np.nan, lw=30, ls='-', c='#77AADD',label="$\mathrm{Rotation}$")
        thrm, = ax.plot(np.nan, np.nan, lw=30, ls='-', c='#CC99BB',label="$\mathrm{Thermohaline}$")
        # semiconv, = ax.plot(np.nan, np.nan, lw=30, ls='-', c='#DDDD77',label="$\mathrm{Semiconvection}$")
        ax.legend(handles=[conv, ov, rot, thrm], loc='lower left', bbox_to_anchor= (0, 1.01), ncol=2,
            borderaxespad=0, frameon=False, handlelength=0.5)


def wrapper(p):
    # set_plotDefaults()
    plt.ioff() # no interactive plots
    fig = plt.figure(figsize=(12,12))
    gs = gridspec.GridSpec(100,100)
    ax = fig.add_subplot(gs[:, :])
    plot_Dmix(p, ax=ax, legend=True)
    ax.set_xlabel(r"$\mathrm{Mass \ } [M_\odot]$")
    ax.set_ylim(10**0, 10**18)
    ax.set_xlim(xmin=0, xmax=21)
    ax.set_ylabel(r"$D \ \mathrm{[cm^2\ s^{-1}]}$")
    ax.axvspan(17, 21, alpha=0.5, color="#808080", zorder=0)
    # ax.set_title(p.split('/')[-1], fontsize=30)
    # plt.tight_layout()
    plt.savefig(root+'/mix_plots/mix_'+f'{profiles.index(p):04.0f}'+'.png')
    plt.close()



if __name__ == "__main__":
    root = "/mnt/home/mrenzo/RUNS/zeta_ophiuchi/M1_25_M2_17_P100_Z0.01/re-run_exponential/"
    LOGS = root+"binary_to_postRLOF/LOGS2/"
    LOGS_single = root+"accretor_to_TAMS/LOGS/"

    profiles_bin = sorted(glob.glob(LOGS+'/profile*.data'), key=get_age_from_profile)
    profiles_single = sorted(glob.glob(LOGS_single+'/profile*.data'), key=get_age_from_profile)
    profiles = profiles_bin+profiles_single

    for p in profiles_bin:
        wrapper(p)
    for p in profiles_single:
        wrapper(p)
    print("done!")
