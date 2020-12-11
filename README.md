# $\zeta$ Ophiuchi

author: [Mathieu Renzo](mailto:mrenzo@flatironinstitute.org)

Project on self-consistent MESA (+CMFGEN?) modeling of [$\zeta$ ophiuchi](http://simbad.u-strasbg.fr/simbad/sim-id?Ident=zeta+ophiuchi&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id).

## $\zeta$ Ophiuchi's observational literature

This doesn't mean to be a comprehensive list:

[Sota et al. 2014](https://ui.adsabs.harvard.edu/abs/2014ApJS..211...10S/abstract),
[Villamariz & Herrero 2005](https://www.aanda.org/articles/aa/pdf/2005/40/aa2848-05.p),
[Howarth & Smith 2001](https://ui.adsabs.harvard.edu/abs/2001MNRAS.327..353H/abstract),
[Tetzlaff al. 2010](https://ui.adsabs.harvard.edu/abs/2010MNRAS.402.2369T/abstract),
[Tetzlaff et al. 2011](https://ui.adsabs.harvard.edu/abs/2011MNRAS.410..190T/abstract)

From its runaway nature we know it is associated to PSR B1929+10 (PSR
J1932+1059) their parent association is probably Upper Scorpio in Sco
OB1 which is likely to be slightly sub-solar Z based on the precision
fitting of pre-MS asteroseismology of one star by [Murphy et al. 2020](https://ui.adsabs.harvard.edu/abs/2020arXiv201111821M/abstract).
Indeed, Z lower than Zsun is required to fit the temperature required
by Villamariz & Herrero 2005.

## Miscellaneous notes

### MESA

Right now MESA evolves a binary until the donor's He core depletion.
If one wants to finish the evolution of the accretor, you should run
it separately as a single star (to avoid the difficult computation of
the very late phases for the donor to be the bottleneck). MESA saves a
model that can be used as initial condition for this, alternatively
one can use a photo. Preliminary tests show this is equivalent.

### Metallicity

The `epsilon(oxygen)` reported in Villamariz & Herrero 2005
corresponds to the Solar value according to Asplund 2005 (cf. Japeli
et al. 2018). There is no way to enhance oxygen through binary
interactions and/or rotational mixing, so this is a good indicator of
the initial Z (or suggests a non-solar mixture)

### Kinematics

Unfortunately the `ruwe` for $\zeta$ Ophiuchi is >4, so EDR3 is
useless. I suspect this is due to the stellar variability and the very
large visual magnitude `Mv = -4.2`.


