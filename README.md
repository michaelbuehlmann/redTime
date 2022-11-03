# redTime

Cosmological perturbation theory code of Upadhye (2019) using the Time-Renormalization Group method of Pietroni (2008) and computing the redshift-space distortion corrections of Taruya, Nishimichi, and Saito (2010)

## Installation (CMake)

redTime requires the [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/). CMake will automatically look
for an installation of GSL. If it cannot be detected automatically, you can set the ``GSL_ROOT_DIR`` environmental
variable (the script expects to find libraries at ``$GSL_ROOT_DIR/lib`` and the GSL headers at
``$GSL_ROOT_DIR/include/gsl``). Check the [CMake documentation](https://cmake.org/cmake/help/latest/module/FindGSL.html)
for more information.



```bash
mkdir build && cd build
cmake ..
make
```

## Running

redTime is run from the command line. The executable can either be directly
called or via a script that reads cosmological parameters from a file, calls
CAMB and then passes the parameters and transferfunctions to the redTime
executable.

### Running via script


### Running redTime directly

The code requires the following files:

* ``params_redTime.dat``, listing the cosmological parameters and code inputs.
  See the example file for the format.
* A CAMB transfer function file (z=0) with 7 columns. The filename is specified
  in ``params_redTime.dat``
* if running with massive neutrinos, then a series of higher-redshift transfer
  functions.

  See ``examples/1_redTime/run.sh`` for an example.

## Outputs

At each redshift requested in params_redTime.dat, redTime produces a 17-column
table of results.  The columns are:

```
1   wave number k [h/Mpc]
2   CDM+Baryon growth factor D normalized to 1 at z=0
3   logarithmic growth derivative d ln(D) / d ln(a)
4   linear CDM+Baryon power spectrum
5   ratio B_nu(k,z) of linear neutrino to CDM+Baryon density
    contrasts, normalized to 1 at z=0
6   logarithmic derivative d ln(B_nu) / d ln(a)
7   linear neutrino power spectrum
8   non-linear CDM+Baryon density-density spectrum
9   non-linear CDM+Baryon density-velocity spectrum, where
    velocity refers to the scalar velocity potential
    theta = -divergence(v) / (a*H)
10  non-linear CDM+Baryon velocity-velocity spectrum
11  mu^2 component of the A(k,mu) RSD correction term from
    Taruya, Nishimichi, and Saito, Phys.Rev. D82 (2010) 063522
    [arXiv:1006.0699] [TNS2010]
12  mu^4 component of A(k,mu) from [TNS2010]
13  mu^6 component of A(k,mu) from [TNS2010]
14  mu^2 component of the B(k,mu) RSD correction from [TNS2010]
15  mu^4 component of B(k,mu) [TNS2010]
16  mu^6 component of B(k,mu) [TNS2010]
17  mu^8 component of B(k,mu) [TNS2010]
```

Preceding each table of results is a comment line beginning with ``#`` and
listing the scale factor and redshift.  The dimensionful quantities ``H``
(Hubble parameter) and ``sigma_v^2`` (linear velocity dispersion length scale)
are also listed, in units of Mpc/h raised to the appropriate power, with the
speed of light c = 1.  (Thus, for example, H at z=0 is always 0.0003336, or
100 km/sec divided by the speed of light.)


## Tuning performance and outputs

Default code performance parameters have been chosen for a combination of speed
and accuracy.   For greater accuracy as well as combatibility with the earlier
redTime version 0.1 (available
[here](http://www.hep.anl.gov/cosmology/pert.html)) at the cost of speed, the
following quantities may be modified prior to compilation of redTime:

  1.  In the function Beta_P(double, input_data) in the file
      AU_cosmological_parameters.h, k_min may be changed to 1e-5 and
      k_max to 20.

  2.  In the function D_dD(double,input_data,double*), also in the file
      AU_cosmological_parameters.h, n_lnk may be changed to 1000 and
      a_early to 1e-50.

  3.  In the file redTime.cc the number of k steps nk may be increased to 256,
      and the size of the extrapolated power spectrum np may be reduced from
      4*nk to 8*nk.

The output format and set of quantities printed may be modified prior to
compilation.  In the file redTime.cc, PREC controls the number of decimal places
printed by the code.  Setting PRINTA, PRINTI, PRINTQ, or PRINTBIAS to 1 includes
among the outputs A_{acd,bef}, I_{acd,bef}, Q^{\ell}_{abc}, or the bias
integrals of McDonald and Roy, JCAP 0908 (2009) 020 [arXiv:0902.0991],
respectively, as summarized in Upadhye, JCAP 1905 (2019) no.05, 041
[arXiv:1707.09354].

In cosmologies with scale-dependent growth factors, such as those with massive
neutrinos, the redshift at which 1-loop mode-coupling integrals are computed may
affect the results.  This redshift is called z1l in the file redTime.cc and is
set to 10 by default.  For compatibility with the earlier redTime version 0.1,
set z1l = C.z_in()) .